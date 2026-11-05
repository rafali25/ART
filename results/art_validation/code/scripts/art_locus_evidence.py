"""Combine family, context, array and completeness evidence; never use spacer length."""


def classify_locus(*, rt_family, partner_family, compatible_context, array_status,
                   copies, upstream_shortfall=0, ambiguous_nt=0, annotation_gap=False):
    reasons = []
    if not rt_family:
        reasons.append('ART_RT_family_not_supported')
    if not partner_family:
        reasons.append('ART_partner_family_not_supported')
    if not compatible_context:
        reasons.append('incompatible_or_unresolved_gene_context')
    if reasons:
        return 'ART_locus_not_supported', reasons
    if upstream_shortfall or ambiguous_nt or annotation_gap:
        if upstream_shortfall: reasons.append('upstream_sequence_incomplete')
        if ambiguous_nt: reasons.append('ambiguous_sequence')
        if annotation_gap: reasons.append('annotation_gap')
        return 'incomplete_ART_locus_assessment', reasons
    if copies == 2:
        return 'exploratory_two_copy_ART_candidate', ['two_copies_only']
    if copies >= 3 and array_status == 'accepted_array_candidate':
        return 'supported_ART_locus_candidate', []
    return 'ART_locus_requires_review', [array_status]
