import sage.all

from .spectral_curves import (
    LocalSpectralCurve,
    RSCurve,
    AiryCurve,
    is_admissible_ramification_profile,
)

from .graph_sums import (
    QAiryGraph,
    zero_qa_graph,
    generate_qa_graphs,
    WGraph,
    zero_W_graph,
    generate_W_graphs,
    generate_W_k_graphs,
    get_poly_dict,
    get_poly_dict_list,
    add_contribution_to_poly_dict,
    add_contribution_of_W_k_graphs,
)

from .generating_function import (
    GeneratingFunction,
    poly_dict_to_str,
    term_to_str,
    load_F_for_rs_curve,
)

from .partitions import (
    partitions_with_min,
    decreasing_partitions_with_min,
    cached_decreasing_red_partitions_with_min,
    decreasing_red_partitions_with_min,
    increasing_partitions_with_min,
    cached_increasing_red_partitions_with_min,
    increasing_red_partitions_with_min,
    tuple_automorphism_number,
    get_pairs,
    subsets_k,
    partition_k,
)

from .misc import (
    r_factorial,
    add_to_dict,
    get_nested_element,
    compute_rs_qairy_index,
    delta,
)
