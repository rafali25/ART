"""Meaningful checks before expanding paid/requested prediction batches."""
import copy
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from email.message import Message
import urllib.error

from art_bindings_common import ROOT,cropped_a3m,msa_csv,polymer,payload_hash
from art_bindings_api import normalize_response,run_request,archive_gateway_failures
import art_bindings_api as transport


class BindingInputs(unittest.TestCase):
    def test_interval_projection_preserves_query_and_homolog_insertions(self):
        self.assertEqual(cropped_a3m('MA--KLT','MAVVK-T',0,5),'MAvvK-T')
        self.assertEqual(cropped_a3m('MA--KLT','MAVVK-T',2,5),'K-T')
        self.assertEqual(cropped_a3m('MA--KLT','MA--KLT',1,4),'AKL')
        with self.assertRaises(ValueError):cropped_a3m('MA','MA',0,3)

    def test_explicit_query_and_unpaired_keys(self):
        query=polymer('A','MAKT','S_query_only')
        self.assertEqual(query['msa']['controlled_input']['a3m']['alignment'],'>query\nMAKT\n')
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'msa.csv';path.write_text(msa_csv(['MAKT','M-KT']))
            unpaired=polymer('A','MAKT','U_requested_independent_MSA_order',path)
            content=unpaired['msa']['controlled_input']['csv']['alignment']
            self.assertTrue(content.startswith('key,sequence\n'))
            self.assertTrue(all(not r['key'] for r in csv.DictReader(io.StringIO(content))))
            path.write_text('sequence,key\nMAKT,shared_locus\n')
            with self.assertRaises(ValueError):polymer('A','MAKT','U_requested_independent_MSA_order',path)

    def test_cache_hash_tracks_sampling_and_alignment_contents(self):
        original=dict(polymers=[polymer('A','MAKT','S_query_only')],diffusion_samples=1)
        multi=dict(original,diffusion_samples=3)
        self.assertNotEqual(payload_hash(original),payload_hash(multi))
        self.assertEqual(payload_hash(original),payload_hash(dict(reversed(list(original.items())))))

    def test_multiple_models_are_sequence_verified_and_duplicates_flagged(self):
        old=ROOT/'data/ART binding/api_pilot_20261105/S_r1_p1'
        body=json.loads((old/'response.json').read_text());payload=json.loads((old/'request.json').read_text())
        payload['diffusion_samples']=3
        for key,value in list(body.items()):
            if isinstance(value,list) and len(value)==1:body[key]=value*3
        with tempfile.TemporaryDirectory() as directory:
            folder=Path(directory)
            result=normalize_response(payload,body,folder)
            self.assertEqual(len(result['models']),3)
            self.assertEqual(result['verified_sequences'],6)
            self.assertFalse(result['structures_have_unique_hashes'])
            broken=copy.deepcopy(payload);broken['polymers'][0]['sequence']='A'+broken['polymers'][0]['sequence'][1:]
            with self.assertRaises(ValueError):normalize_response(broken,body,folder)


class TransportRecovery(unittest.TestCase):
    def test_accepted_request_is_polled_without_duplicate_post_after_gateway_error(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);payload=dict(diffusion_samples=3,polymers=[])
            digest=payload_hash(payload);request=root/'request.json';request.write_text(json.dumps(payload))
            record=dict(request_hash=digest,request_path='request.json',stages=['test'])
            class Response:
                def __init__(self,status,body,request_id=None):
                    self.status=status;self.body=body;self.headers=Message()
                    if request_id:self.headers['NVCF-REQID']=request_id
                def __enter__(self):return self
                def __exit__(self,*args):pass
                def read(self,*args):
                    # The accepted invocation must already be recoverable even
                    # if reading/decoding the body fails immediately afterwards.
                    state=json.loads((root/'runs'/digest/'state.json').read_text())
                    if self.headers.get('NVCF-REQID'):assert state['request_id']=='request-test'
                    return self.body
            headers=Message();headers['NVCF-REQID']='request-test'
            gateway=urllib.error.HTTPError('https://health.api.nvidia.com/v1/status/request-test',504,'Gateway',headers,io.BytesIO(b''))
            calls=[];responses=iter([Response(202,b'{}','request-test'),gateway,Response(200,b'{}','request-test')])
            class Opener:
                def open(self,request,**kwargs):
                    calls.append(request)
                    item=next(responses)
                    if isinstance(item,Exception):raise item
                    return item
            with patch.object(transport,'OUT',root),patch.object(transport.urllib.request,'build_opener',return_value=Opener()),patch.object(transport.time,'sleep'),patch.object(transport,'normalize_response',return_value={'models':[{}, {}, {}]}):
                state=run_request(record,'test-credential',{'endpoint':'https://health.api.nvidia.com/test'})
            self.assertEqual(state['status'],'complete')
            self.assertEqual([c.method for c in calls],['POST','GET','GET'])
            self.assertEqual(calls[0].get_header('Nvcf-poll-seconds'),'300')
            self.assertEqual(json.loads((root/'runs'/digest/'http_error_001.json').read_text())['headers']['NVCF-REQID'],'request-test')

    def test_gateway_retry_is_archived_but_unknown_socket_submission_is_not_retried(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);records=[]
            for digest,status in [('gateway','failed'),('socket','submission_uncertain')]:
                folder=root/'runs'/digest;folder.mkdir(parents=True)
                (folder/'state.json').write_text(json.dumps(dict(status=status,request_id=None)))
                (folder/'http_error_000.json').write_text(json.dumps(dict(status=504)))
                records.append(dict(request_hash=digest))
            with patch.object(transport,'OUT',root):revisions=archive_gateway_failures(records)
            self.assertEqual(len(revisions),1)
            self.assertTrue((root/'runs/gateway/attempts/001/state.json').exists())
            self.assertFalse((root/'runs/gateway/state.json').exists())
            self.assertTrue((root/'runs/socket/state.json').exists())


if __name__=='__main__':unittest.main()
