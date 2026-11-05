import unittest
import json
import numpy as np
from art_binding_msa_control import separate_chains,audit_msa


class InputControls(unittest.TestCase):
    def test_split_preserves_each_homolog_and_deletion(self):
        msa=np.array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]])
        deletion=np.arange(12).reshape(4,3);paired=np.ones_like(msa);chains=np.array([0,0,1,1])
        m,d,p=separate_chains(msa,deletion,paired,chains,0)
        self.assertEqual(m.shape,(4,5));np.testing.assert_array_equal(m[:,0],msa[:,0])
        for c in [0,1]:
            cols=[i for i in range(1,m.shape[1]) if np.any(m[chains==c,i])]
            np.testing.assert_array_equal(m[chains==c][:,cols],msa[chains==c,1:])
            np.testing.assert_array_equal(d[chains==c][:,cols],deletion[chains==c,1:])
        result=audit_msa(m,p,chains,0,'U')
        self.assertEqual(result['nonquery_multiple_chain_rows'],0)
        self.assertEqual(json.loads(json.dumps(result))['rows'],5)

    def test_single_query_is_retained(self):
        msa=np.array([[1],[2],[3],[4]]);chains=np.array([0,0,1,1])
        m,d,p=separate_chains(msa,np.zeros_like(msa),np.ones_like(msa),chains,0)
        np.testing.assert_array_equal(m,msa)
        self.assertEqual(audit_msa(m,p,chains,0,'S')['rows'],1)

    def test_audit_rejects_hidden_pairing_and_homologs_in_S(self):
        msa=np.array([[1,2],[3,4]]);chains=np.array([0,1])
        with self.assertRaises(ValueError):audit_msa(msa,np.zeros_like(msa),chains,0,'U')
        m,d,p=separate_chains(msa,np.zeros_like(msa),np.ones_like(msa),chains,0)
        with self.assertRaises(ValueError):audit_msa(m,p,chains,0,'S')

    def test_incompatible_shapes_are_rejected(self):
        with self.assertRaises(ValueError):separate_chains(np.zeros((2,2)),np.zeros((2,3)),np.zeros((2,2)),[0,1],0)


if __name__=='__main__':unittest.main()
