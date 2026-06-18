//! Cross-language parity check: the Rust `neighbours` must match the JS demo
//! (web/codec.js) exactly, so the two implementations interoperate.
use rvdna::storage::fountain::neighbours;

#[test]
fn neighbours_match_js_reference() {
    // (seed, k) => (degree, idxs) taken from web/codec.js neighbours().
    let cases: &[((u32, usize), (u32, &[usize]))] = &[
        ((0, 1), (1, &[0])),
        ((0, 2), (2, &[0, 1])),
        ((12648430, 16), (6, &[2, 4, 7, 9, 12, 15])),
        ((12648454, 16), (16, &[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])),
        ((7, 20), (2, &[4, 17])),
        ((99, 50), (2, &[8, 14])),
        ((123456, 19), (2, &[7, 10])),
    ];
    for ((seed, k), (deg, idxs)) in cases {
        let (d, got) = neighbours(*seed, *k);
        assert_eq!(d, *deg, "degree mismatch seed={seed} k={k}");
        assert_eq!(got, idxs.to_vec(), "idxs mismatch seed={seed} k={k}");
    }
}
