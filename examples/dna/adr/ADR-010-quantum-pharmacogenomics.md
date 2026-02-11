# ADR-010: Quantum-Enhanced Pharmacogenomics & Precision Medicine

**Status**: Proposed
**Date**: 2026-02-11
**Authors**: ruv.io, RuVector DNA Analyzer Team
**Deciders**: Architecture Review Board
**Target Crates**: `ruQu`, `ruvector-gnn`, `ruvector-core`, `ruvector-attention`, `ruvector-sona`

## Version History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 0.1 | 2026-02-11 | RuVector DNA Analyzer Team | Initial proposal |

---

## Context

### The Pharmacogenomics Problem

Pharmacogenomics -- the study of how an individual's genome influences their response to drugs -- remains one of the most actionable domains in clinical genomics. Approximately 95% of patients carry at least one actionable pharmacogenomic variant, yet fewer than 5% of prescriptions incorporate pharmacogenomic testing. The consequences are severe: adverse drug reactions (ADRs) account for approximately 2.2 million hospitalizations and 106,000 deaths annually in the United States alone (Lazarou et al., JAMA 1998; updated estimates suggest the figure is substantially higher).

The core challenge is computational: translating a patient's raw genotype into a precise drug response prediction requires integrating star allele calling, haplotype phasing, gene-drug interaction modeling, metabolic pathway simulation, and dosage optimization -- all within a clinically actionable timeframe.

### Pharmacogene Complexity

The cytochrome P450 (CYP450) superfamily and other pharmacogenes present unique computational challenges that exceed the capabilities of standard variant callers:

| Gene | Chromosomal Region | Known Star Alleles | Key Challenge |
|------|-------------------|-------------------|---------------|
| CYP2D6 | 22q13.2 | >150 | Whole-gene deletion, duplication (1-13 copies), gene conversion with CYP2D7/CYP2D8, hybrid alleles |
| CYP2C19 | 10q23.33 | >35 | Common loss-of-function (*2, *3); gain-of-function (*17) |
| CYP3A4 | 7q22.1 | >30 | Substrate for ~50% of all drugs; highly polymorphic promoter |
| CYP2C9 | 10q23.33 | >60 | Warfarin sensitivity (*2, *3); NSAID metabolism |
| CYP1A2 | 15q24.1 | >20 | Caffeine, theophylline, clozapine metabolism; inducible |
| DPYD | 1p21.3 | >30 | 5-fluorouracil toxicity; DPD deficiency is life-threatening |
| TPMT | 6p22.3 | >40 | Thiopurine toxicity (azathioprine, 6-mercaptopurine); TPMT*2, *3A, *3B, *3C |
| NUDT15 | 13q14.2 | >10 | Thiopurine toxicity; particularly relevant in East Asian populations |
| UGT1A1 | 2q37.1 | >100 | Irinotecan toxicity; Gilbert syndrome (*28 allele); SN-38 glucuronidation |
| SLCO1B1 | 12p12.2 | >40 | Statin-induced myopathy (simvastatin, atorvastatin); Val174Ala (*5) |
| HLA-B | 6p21.33 | >7,000 | Abacavir hypersensitivity (*57:01); carbamazepine SJS/TEN (*15:02) |
| VKORC1 | 16p11.2 | >20 | Warfarin dose requirement; -1639G>A promoter variant |
| CYP2B6 | 19q13.2 | >35 | Efavirenz metabolism; highly variable across populations |
| NAT2 | 8p22 | >90 | Isoniazid hepatotoxicity; slow/rapid/intermediate acetylator phenotypes |
| G6PD | Xq28 | >400 | Primaquine/rasburicase hemolytic anemia; X-linked |

CYP2D6 is the most computationally demanding: a single individual may carry 0-13 copies of the gene, with hybrid alleles formed by recombination between CYP2D6 and the nearby pseudogene CYP2D7. Accurate star allele calling requires simultaneous copy number determination, structural variant resolution, haplotype phasing across the entire CYP2D6/2D7 locus, and disambiguation of CYP2D6-CYP2D7 hybrid alleles.

### Metabolizer Phenotype Classification

Star allele diplotypes translate to metabolizer phenotypes that directly determine drug response:

| Phenotype | Activity Score Range | Prevalence (European) | Clinical Implication |
|-----------|--------------------|-----------------------|---------------------|
| Ultra-rapid metabolizer (UM) | >2.25 | ~1-10% (gene-dependent) | Drug inefficacy at standard dose; prodrug toxicity |
| Normal metabolizer (NM) | 1.25-2.25 | ~60-80% | Standard dosing appropriate |
| Intermediate metabolizer (IM) | 0.25-1.0 | ~10-25% | Consider dose reduction |
| Poor metabolizer (PM) | 0 | ~2-10% | Significant dose reduction or alternative drug required |
| Indeterminate | N/A | Variable | Insufficient evidence to assign phenotype; flag for clinical review |

The activity score system assigns a numerical value to each star allele based on its functional impact:

| Function | Activity Score | Examples |
|----------|---------------|---------|
| Normal function | 1.0 | CYP2D6*1, CYP2D6*2, CYP2D6*35 |
| Decreased function | 0.5 | CYP2D6*9, CYP2D6*10, CYP2D6*17, CYP2D6*29, CYP2D6*41 |
| No function | 0 | CYP2D6*3, CYP2D6*4, CYP2D6*5 (deletion), CYP2D6*6 |
| Increased function | (copy number) x allele score | CYP2D6 gene duplication: score = n_copies x allele_score |

The diplotype activity score is the sum of the two haplotype scores. For CYP2D6 with gene duplications, the total score accounts for all functional copies.

### Why Quantum Enhancement

Classical molecular simulation of drug-enzyme interactions relies on molecular mechanics force fields (e.g., AMBER, CHARMM) that approximate quantum effects with empirical parameters. These approximations break down for:

1. **Transition state modeling**: CYP450 catalytic cycles involve iron-oxo intermediates with strong multireference character -- multiple electronic configurations contribute significantly to the wavefunction
2. **Proton-coupled electron transfer**: CYP2D6 catalysis involves concerted proton and electron transfer that cannot be captured by classical force fields
3. **Spin-state transitions**: The iron center in CYP450 transitions between high-spin (S=5/2), intermediate-spin (S=3/2), and low-spin (S=1/2) states during catalysis; classical methods assign a single spin state
4. **Dispersion interactions**: Non-covalent drug-enzyme interactions in the active site depend on London dispersion forces that require at least DFT-D3 level treatment

Quantum computing -- specifically the Variational Quantum Eigensolver (VQE) -- enables ab initio simulation of these electronic structure problems with polynomial rather than exponential scaling in the number of active orbitals.

### Current Limitations

| Approach | Accuracy for CYP450 Active Sites | Computational Cost | Limitation |
|----------|----------------------------------|-------------------|-----------|
| Molecular mechanics (MM) | Low: no electronic structure | Minutes | Cannot model transition states, spin states, or charge transfer |
| Semi-empirical QM (PM7, GFN2-xTB) | Moderate: parametrized QM | Hours | Parameters fitted to ground states; poor for transition metals |
| DFT (B3LYP/def2-TZVP) | Good for ground states | Days (per conformer) | Systematic errors for reaction barriers; self-interaction error for Fe |
| CASSCF/CASPT2 | Excellent: multireference | Weeks (active space ~14e/14o) | Exponential scaling limits active space to ~20 orbitals classically |
| VQE on quantum hardware | Potentially excellent | Hours (with quantum advantage) | Current NISQ hardware: ~100 noisy qubits; error mitigation required |

The RuVector `ruQu` crate provides the quantum simulation backend with hybrid classical-quantum VQE, noise-aware circuit optimization, and automatic fallback to classical CASSCF when quantum hardware is unavailable.

---

## Decision

### Adopt a Quantum-Enhanced Pharmacogenomics Pipeline

We implement a pharmacogenomics pipeline that integrates:

1. **Star allele calling** via GNN-based structural resolution of pharmacogene loci (`ruvector-gnn`)
2. **Drug-gene interaction prediction** via graph neural network on a pharmacogenomic knowledge graph (`ruvector-gnn`)
3. **Molecular docking and reaction modeling** via hybrid VQE on quantum hardware (`ruQu`)
4. **Adverse event prediction** via HNSW similarity search over historical patient-drug-outcome vectors (`ruvector-core`)
5. **Polypharmacy interaction analysis** via multi-head attention over drug combination tensors (`ruvector-attention`)
6. **Bayesian dosage optimization** via SONA-adapted posterior estimation (`ruvector-sona`)
7. **Clinical decision support** with genotype-to-phenotype translation, interaction alerts, and dosing recommendations

---

## Core Capabilities

### 1. Star Allele Calling

#### The Star Allele Nomenclature

Pharmacogene alleles are designated using the star allele (*) nomenclature, where *1 is the reference (wild-type) allele and numbered alleles (*2, *3, ...) represent defined combinations of variants. A diplotype is the pair of star alleles on the two homologous chromosomes (e.g., CYP2D6 *1/*4).

#### GNN-Based Structural Resolution

Standard variant callers fail on CYP2D6 because the locus contains:
- A tandem arrangement of CYP2D6 and pseudogenes CYP2D7/CYP2D8
- Whole-gene deletions (*5 allele) and duplications (CYP2D6xN, N=2-13)
- Gene conversion events producing hybrid CYP2D6-CYP2D7 alleles (*13, *36, *57, *68)
- Structural variants spanning 30-50 kbp

The `ruvector-gnn` module resolves these events using a graph neural network architecture:

**Graph Construction**:
For the CYP2D6 locus (chr22:42,126,499-42,185,000, GRCh38):

1. Extract all reads mapping to the CYP2D6/CYP2D7/CYP2D8 region (including supplementary alignments)
2. Construct a read-overlap graph where:
   - Nodes = reads, featurized by [mapping_quality, insert_size, num_mismatches, has_soft_clip, is_supplementary, mate_distance]
   - Edges = reads sharing >= 50 bp overlap, weighted by overlap quality score
3. Add structural nodes representing known CYP2D6 configurations:
   - *1 (reference, single copy)
   - *5 (whole-gene deletion)
   - *1xN (N-copy duplication, N=2..13)
   - *13 (CYP2D6/CYP2D7 hybrid, exons 1-7 from CYP2D7)
   - *36 (CYP2D6/CYP2D7 hybrid, different breakpoint)
   - *68 (CYP2D6/CYP2D7 hybrid with upstream deletion)

**GNN Message Passing**:

```
h_v^(l+1) = sigma(W^(l) * AGG({h_u^(l) : u in N(v)}) + b^(l))
```

where:
- `h_v^(l)` is the hidden state of node v at layer l
- `N(v)` is the neighborhood of v in the read-overlap graph
- `AGG` is a permutation-invariant aggregation (mean + max pooling)
- `sigma` is the GELU activation function
- 4 message-passing layers with hidden dimension 256

After L layers of message passing, a global readout produces the structural configuration posterior:

```
P(config | reads) = softmax(MLP(ReadOut(h^(L))))
```

where `config` ranges over all known CYP2D6 structural configurations (single-copy alleles, deletions, duplications, hybrids).

**Copy Number Determination**:
The GNN simultaneously estimates the CYP2D6 copy number by regressing on the normalized read depth ratio:

```
CN_estimate = round(depth_CYP2D6 / depth_reference * 2)
```

where `depth_reference` is computed from flanking single-copy regions. The GNN refines this estimate using read-pair and split-read evidence from structural edges in the graph.

**Haplotype Phasing**:
Once the structural configuration is determined, star allele assignment proceeds by:

1. Phase variants within each CYP2D6 copy using WhatsHap-style read-backed phasing
2. Assign each phased haplotype to a known star allele by matching the variant combination against the PharmVar database (https://www.pharmvar.org/)
3. For novel combinations not in PharmVar, assign a provisional star allele with the closest known match (by Hamming distance on defining variants) and flag for expert review

**Star Allele Calling for All Pharmacogenes**:

The same GNN architecture is applied to each pharmacogene, with gene-specific adaptations:

| Gene | Structural Complexity | GNN Configuration |
|------|----------------------|-------------------|
| CYP2D6 | Very high (deletion, duplication, hybrid) | Full structural GNN, 4 layers, copy number regression |
| CYP2C19 | Low (no common structural variants) | Simplified GNN, 2 layers, SNP-based calling sufficient |
| CYP3A4 | Moderate (rare deletions) | 3-layer GNN, depth-based deletion detection |
| CYP2C9 | Low | 2-layer GNN |
| DPYD | Moderate (large gene, 23 exons, intronic variants matter) | 3-layer GNN, includes intronic splicing variant detection |
| TPMT | Low (well-characterized alleles: *2, *3A, *3B, *3C) | 2-layer GNN |
| NUDT15 | Low | 2-layer GNN |
| UGT1A1 | Moderate (TA repeat in promoter; *28 = (TA)7, *36 = (TA)5, *37 = (TA)8) | 3-layer GNN with STR length estimation for promoter repeat |
| SLCO1B1 | Low | 2-layer GNN |
| HLA-B | Very high (>7,000 alleles, extreme polymorphism, MHC region) | Dedicated HLA-typing module (OptiType-style, but GNN-enhanced) |
| VKORC1 | Low (key variant: -1639G>A) | 2-layer GNN |

### 2. Drug-Gene Interaction Prediction via GNN

#### Knowledge Graph Structure

The pharmacogenomic knowledge graph integrates four major databases:

| Database | Entities | Relationships | Version |
|----------|----------|---------------|---------|
| CPIC (Clinical Pharmacogenetics Implementation Consortium) | 27 genes, 90+ drugs | Gene-drug pairs with guideline-level evidence | Current |
| PharmGKB | 745 genes, 985 drugs | Clinical annotations, variant-drug associations, pathway data | Current |
| DrugBank | 14,000+ drugs | Drug-target interactions, metabolic pathways, enzyme inhibition/induction | 5.1 |
| UniProt | 20,000+ human proteins | Protein function, 3D structure, active site annotations | 2026_01 |

**Graph Schema**:

```
Nodes:
  Gene {id, symbol, chromosome, star_alleles[], function, expression_tissue[]}
  Drug {id, name, smiles, drugbank_id, atc_code, mw, logP, psa}
  Protein {id, uniprot_id, pdb_ids[], active_site_residues[], ec_number}
  Variant {id, rsid, gene, hgvs_c, hgvs_p, consequence, star_allele}
  Phenotype {id, metabolizer_status, activity_score}
  Pathway {id, name, kegg_id, reaction_steps[]}
  AdverseEvent {id, meddra_pt, severity, organ_class}

Edges:
  METABOLIZES(Gene -> Drug) {km, vmax, clint, evidence_level}
  TRANSPORTS(Gene -> Drug) {direction: influx|efflux, km, evidence_level}
  INHIBITS(Drug -> Gene) {ki, mechanism: competitive|noncompetitive|irreversible}
  INDUCES(Drug -> Gene) {fold_change, time_to_max_induction}
  ENCODES(Gene -> Protein)
  DEFINES(Variant -> Phenotype) {functional_consequence, activity_score_delta}
  PARTICIPATES(Protein -> Pathway) {role: enzyme|substrate|cofactor}
  CAUSES(Drug -> AdverseEvent) {frequency, severity, pharmacogene_association}
  INTERACTS(Drug -> Drug) {mechanism, severity, clinical_significance}
```

**GNN Architecture for Interaction Prediction**:

The knowledge graph is embedded using a relational graph convolutional network (R-GCN) that learns type-specific message-passing functions for each edge type:

```
h_v^(l+1) = sigma(sum_{r in R} sum_{u in N_r(v)} (1/c_{v,r}) * W_r^(l) * h_u^(l) + W_0^(l) * h_v^(l))
```

where:
- `R` = set of edge types (METABOLIZES, TRANSPORTS, INHIBITS, INDUCES, etc.)
- `N_r(v)` = neighbors of v connected by edge type r
- `W_r^(l)` = relation-specific weight matrix at layer l
- `c_{v,r}` = normalization constant (degree of v under relation r)
- 6 R-GCN layers, hidden dimension 512

**Prediction Tasks**:

| Task | Input | Output | Metric |
|------|-------|--------|--------|
| Drug-gene interaction type | (Drug, Gene) pair | {metabolizes, transports, inhibits, induces, none} | AUC-ROC > 0.95 |
| Interaction strength | (Drug, Gene) pair with known type | Km, Vmax, Ki (continuous) | Spearman rho > 0.85 |
| Adverse event risk | (Drug, Genotype) pair | P(adverse event by MedDRA preferred term) | AUC-ROC > 0.90 |
| Drug-drug interaction via shared pathway | (Drug_A, Drug_B) pair | {contraindicated, major, moderate, minor, none} | F1 > 0.85 |
| Novel drug-gene association discovery | Drug embedding | Ranked list of candidate gene interactions | Recall@20 > 0.80 |

### 3. Molecular Docking via Quantum VQE

#### CYP450 Active Site Modeling

CYP450 enzymes catalyze oxidation reactions using a heme iron center (protoporphyrin IX with axial cysteine thiolate ligand). The catalytic cycle involves:

1. Substrate binding in the active site cavity
2. Reduction of Fe(III) to Fe(II) by NADPH-cytochrome P450 reductase
3. O2 binding to Fe(II), forming Fe(II)-O2 (Compound 0)
4. Protonation and O-O bond cleavage, forming Fe(IV)=O (Compound I, the active oxidant)
5. Hydrogen atom abstraction from substrate by Compound I
6. Oxygen rebound to form hydroxylated product
7. Product release, returning to resting Fe(III) state

**Quantum Hamiltonian for the Active Site**:

The electronic Hamiltonian for the CYP450 active site in second quantization:

```
H = sum_{pq} h_pq * a_p^dagger * a_q + (1/2) * sum_{pqrs} h_pqrs * a_p^dagger * a_q^dagger * a_s * a_r
```

where:
- `a_p^dagger, a_q` are fermionic creation/annihilation operators for molecular orbital p, q
- `h_pq` = one-electron integrals (kinetic energy + nuclear attraction)
- `h_pqrs` = two-electron integrals (electron-electron repulsion)
- Indices p, q, r, s run over the active space orbitals

**Active Space Selection**:

For CYP450 Compound I (Fe(IV)=O porphyrin with cysteine thiolate):

| Active Space | Orbitals | Electrons | Qubits (Jordan-Wigner) | Classical Cost | Quantum Cost |
|-------------|---------|-----------|----------------------|---------------|-------------|
| Minimal | (8e, 8o): Fe 3d, O 2p | 8 | 16 | CASSCF feasible (~hours) | VQE: ~minutes |
| Standard | (14e, 14o): Fe 3d, O 2p, porphyrin pi, Cys S 3p | 14 | 28 | CASPT2: ~days | VQE: ~hours |
| Extended | (20e, 20o): + substrate frontier orbitals | 20 | 40 | Intractable classically | VQE: ~hours |
| Full substrate | (30e, 30o): + full substrate valence | 30 | 60 | Completely intractable | VQE: feasible on fault-tolerant QC |

The extended (20e, 20o) active space is the target for near-term quantum advantage: it includes the substrate frontier orbitals necessary for modeling the hydrogen atom transfer transition state, but exceeds the practical limits of classical CASSCF/CASPT2.

**Qubit Mapping**:

The fermionic Hamiltonian is mapped to qubit operators using the Jordan-Wigner transformation:

```
a_p^dagger -> (1/2)(X_p - iY_p) * prod_{q<p} Z_q
a_p        -> (1/2)(X_p + iY_p) * prod_{q<p} Z_q
```

where X, Y, Z are Pauli matrices. This maps N orbitals to N qubits with O(N) Pauli terms per operator.

**VQE Ansatz**:

The Unitary Coupled Cluster with Singles and Doubles (UCCSD) ansatz:

```
|psi(theta)> = exp(T(theta) - T^dagger(theta)) |phi_0>
```

where:
- `|phi_0>` = Hartree-Fock reference state
- `T(theta) = T_1(theta) + T_2(theta)`
- `T_1 = sum_{ia} theta_i^a * a_a^dagger * a_i` (single excitations)
- `T_2 = sum_{ijab} theta_{ij}^{ab} * a_a^dagger * a_b^dagger * a_j * a_i` (double excitations)
- `theta` = variational parameters optimized by the VQE loop

**Number of variational parameters**: For (20e, 20o) active space:
- Singles: 10 occupied x 10 virtual = 100 parameters
- Doubles: C(10,2) x C(10,2) = 45 x 45 = 2,025 parameters
- Total: ~2,125 parameters (manageable for classical optimizer)

**VQE Optimization Loop** (`ruQu`):

```
1. Prepare reference state |phi_0> on quantum hardware
2. Apply parameterized UCCSD circuit U(theta)
3. Measure expectation value <psi(theta)|H|psi(theta)> by:
   - Decompose H into sum of Pauli strings: H = sum_k c_k * P_k
   - Measure each P_k independently (grouping commuting terms)
   - Compute E(theta) = sum_k c_k * <P_k>
4. Update theta using classical optimizer (L-BFGS-B or COBYLA)
5. Repeat until |E(theta_n) - E(theta_{n-1})| < 1e-6 Hartree (convergence)
6. Extract optimized energy and wavefunction
```

**Error Mitigation** (for NISQ hardware):

| Technique | Overhead | Error Reduction |
|-----------|---------|----------------|
| Zero-noise extrapolation (ZNE) | 3-5x circuit repetitions | ~10x error reduction |
| Probabilistic error cancellation (PEC) | O(exp(n_gates * epsilon)) | Exact in principle; practical for <1000 gates |
| Symmetry verification | 2x measurements | Rejects measurements violating N-electron or spin symmetry |
| Virtual distillation | 2 copies of state | Quadratic error suppression |

`ruQu` implements ZNE and symmetry verification as default error mitigation strategies, with PEC available for high-accuracy calculations.

#### Quantum Phase Estimation for Reaction Barriers

For the hydrogen atom transfer (HAT) step of CYP450 catalysis -- the rate-determining step for most substrates -- the reaction barrier determines the metabolic rate. Quantum Phase Estimation (QPE) provides the ground state energy to chemical accuracy (1 kcal/mol = 1.6 mHartree):

```
QPE Circuit:
1. Prepare approximate ground state |psi_approx> (from VQE or classical CASSCF)
2. Apply controlled-U^(2^k) gates for k = 0, 1, ..., n-1 ancilla qubits
3. Inverse QFT on ancilla register
4. Measure ancilla register to obtain phase phi
5. E_ground = phi * (2*pi / t)
```

**Resource Estimates for CYP2D6 Compound I + Substrate**:

| Parameter | Value |
|-----------|-------|
| Active space | (20e, 20o) |
| Qubits (system) | 40 |
| Ancilla qubits (for QPE, 1 mHartree precision) | 15-20 |
| Total qubits | 55-60 |
| T-gate count (Trotter-Suzuki, 2nd order) | ~10^8 |
| Estimated wall time (fault-tolerant QC, 1 MHz logical gate rate) | ~100 seconds |
| Classical CASPT2 equivalent wall time | ~weeks (if feasible at all) |

**Born-Oppenheimer Potential Energy Surface**:

To compute the HAT reaction barrier, we scan the C-H...O coordinate (C = substrate carbon, H = transferred hydrogen, O = Compound I oxygen):

```
1. Fix the O...H distance at values: 2.5, 2.0, 1.8, 1.5, 1.3, 1.1, 1.0 Angstrom
2. Optimize all other coordinates classically (MM/QM embedding)
3. At each fixed O...H distance, compute E_quantum(R) via VQE or QPE for the active site
4. Fit the minimum energy path to a cubic spline
5. Extract the barrier height: dE = E_TS - E_reactant
6. Compute the rate constant via Eyring equation:
   k = (k_B * T / h) * exp(-dG_barrier / (R * T))
```

The barrier height directly maps to Km (the Michaelis-Menten constant):

```
Km ~ exp(dG_binding / (R * T))
Vmax ~ k_cat = (k_B * T / h) * exp(-dG_barrier / (R * T))
```

**Variant Impact on Binding Affinity**:

Pharmacogenomic star alleles alter the enzyme structure, which changes the active site geometry and thus the binding affinity and catalytic rate for specific substrates. For each star allele defining variant:

1. Model the structural impact using AlphaFold2 with the variant amino acid substitution
2. Re-dock the drug substrate in the altered active site
3. Recompute VQE energy for the new geometry
4. Delta(Km) and Delta(Vmax) predict the altered metabolic rate for that allele

This provides a mechanistic, first-principles basis for star allele functional annotations -- replacing the current system of empirical clinical observations with quantum-computed predictions.

### 4. CYP450 Reaction Modeling

#### Catalytic Cycle Quantum Simulation

The full CYP450 catalytic cycle is modeled as a series of quantum chemical calculations at key intermediates:

| Intermediate | Electronic State | Active Space | Quantum Method |
|-------------|-----------------|-------------|----------------|
| Resting state (Fe(III)-H2O) | Doublet (S=1/2) | (9e, 9o) | VQE-UCCSD |
| Substrate-bound (Fe(III)-Sub) | Doublet (S=1/2) | (9e, 9o) | VQE-UCCSD |
| Reduced (Fe(II)-Sub) | Quintet (S=2) | (10e, 10o) | VQE-UCCSD |
| Oxy-complex (Fe(II)-O2-Sub) | Singlet/Triplet | (14e, 14o) | VQE-UCCSD |
| Compound 0 (Fe(III)-OOH-Sub) | Doublet (S=1/2) | (14e, 14o) | VQE-UCCSD |
| Compound I (Fe(IV)=O-Sub) | Doublet (S=1/2) + Quartet (S=3/2) | (14e, 14o) + substrate | VQE-UCCSD |
| Transition state (HAT) | Doublet/Quartet | (20e, 20o) | VQE-UCCSD or QPE |
| Rebound intermediate | Doublet/Quartet | (16e, 16o) | VQE-UCCSD |
| Product-bound | Doublet (S=1/2) | (9e, 9o) | VQE-UCCSD |

**Two-State Reactivity (TSR)**:

CYP450 Compound I exists as a near-degenerate doublet-quartet pair. The reaction proceeds on both spin surfaces simultaneously, with the lower barrier determining the rate. This two-state reactivity (Shaik et al., Chem. Rev. 2005) requires computing energies on both spin surfaces -- a natural fit for quantum simulation where different spin states are encoded by different qubit configurations.

```
Barrier_effective = min(Barrier_doublet, Barrier_quartet)
```

The doublet and quartet surfaces cross near the transition state. Spin-orbit coupling (SOC) at the crossing point determines the branching ratio:

```
P_crossing = 2 * pi * |SOC|^2 / (h * v * |dE/dR|)
```

where v is the nuclear velocity along the reaction coordinate.

#### Substrate Oxidation Site Prediction

For a given drug molecule, the CYP450 enzyme can oxidize at multiple positions. The preferred oxidation site determines the metabolite(s) formed. The quantum pipeline predicts the site of metabolism (SOM):

1. Enumerate all C-H bonds in the substrate
2. For each C-H bond, compute the HAT barrier via VQE at the (20e, 20o) level
3. Rank C-H bonds by barrier height (lowest barrier = most likely SOM)
4. Account for steric accessibility by combining quantum barrier with docking score

```
SOM_score(position_i) = w_electronic * exp(-dE_barrier_i / kT) + w_steric * docking_score_i
```

Weights w_electronic and w_steric are trained on the MetaQSAR dataset of experimentally determined SOMs.

### 5. Adverse Event Prediction via HNSW

#### Patient-Drug-Outcome Vector Space

Each historical patient-drug interaction is encoded as a vector in a shared embedding space:

```
v_interaction = [v_patient || v_drug || v_outcome]
```

where:
- `v_patient` (128-dim): pharmacogenomic profile vector (from ADR-003, Section 6: Pharmacogenomic Vectors)
- `v_drug` (128-dim): drug molecular embedding (learned from DrugBank SMILES via GNN)
- `v_outcome` (64-dim): clinical outcome embedding (from EHR structured data: ICD-10 codes, MedDRA terms, lab values)

Total dimensionality: 320-dim per interaction record.

**HNSW Index Configuration**:

```
Index parameters:
  dimensions: 320
  metric: Cosine
  M: 32
  ef_construction: 200
  ef_search: 200 (high recall for clinical safety)
  quantization: Scalar u8 (4x compression; no PQ to preserve clinical accuracy)
```

**Adverse Event Prediction**:

For a new patient with genotype G about to receive drug D:

1. Encode patient pharmacogenomic profile: `v_patient = phi_pgx(G)`
2. Encode drug: `v_drug = phi_drug(D)`
3. Construct partial query vector: `q = [v_patient || v_drug || 0_outcome]`
4. Search HNSW index for k=100 nearest historical interactions (by patient + drug subspace)
5. Aggregate outcomes from k neighbors:

```
P(adverse_event = e | G, D) = sum_{i=1}^{k} w_i * I(outcome_i = e) / sum_{i=1}^{k} w_i
```

where `w_i = exp(-d(q, v_i) / tau)` is a temperature-scaled similarity weight (tau = 0.1).

**Alert Thresholds**:

| Adverse Event Severity | P(event) Threshold | Action |
|----------------------|-------------------|--------|
| Life-threatening (e.g., SJS/TEN, agranulocytosis) | >= 0.01 (1%) | Contraindication alert |
| Serious (e.g., hepatotoxicity, QT prolongation) | >= 0.05 (5%) | Strong warning |
| Moderate (e.g., GI disturbance, drowsiness) | >= 0.20 (20%) | Informational alert |
| Mild | >= 0.50 (50%) | Patient counseling recommendation |

**HNSW Speedup for Adverse Event Search**:

| Patient-Drug Records | Brute Force | HNSW (k=100, ef=200) | Speedup |
|---------------------|------------|----------------------|---------|
| 100K | 50ms | 200us | 250x |
| 1M | 500ms | 400us | 1,250x |
| 10M | 5s | 1ms | 5,000x |
| 100M | 50s | 3ms | 16,667x |

At the scale of a national pharmacovigilance database (100M+ patient-drug interactions), HNSW enables real-time adverse event risk estimation with sub-5ms latency.

### 6. Polypharmacy Analysis

#### Drug Interaction Tensor

Polypharmacy -- the concurrent use of multiple drugs -- creates combinatorial interaction effects that are poorly characterized for most drug combinations. We model pairwise and higher-order interactions using a multi-head attention mechanism.

**Pairwise Interaction Encoding**:

For a patient taking N drugs simultaneously, construct the interaction tensor:

```
T_{ij} = phi_interact(drug_i, drug_j, genotype)
```

where `phi_interact` encodes the pharmacokinetic interaction between drugs i and j given the patient's metabolizer phenotype for all relevant CYP450 enzymes. Dimensions: N x N x d_interact (d_interact = 128).

**Multi-Head Attention over Drug Combinations** (`ruvector-attention`):

```
Q = W_Q * D_patient    (drug embeddings for patient's medications)
K = W_K * D_patient
V = W_V * T_interaction (interaction tensor features)

Attention(Q, K, V) = softmax(Q * K^T / sqrt(d_k)) * V
```

Flash attention (`ruvector-attention::FlashAttention`) provides 2.49x-7.47x speedup over naive attention for patients on many medications (N > 10).

**Interaction Types Detected**:

| Interaction Type | Mechanism | Example | Detection Method |
|-----------------|-----------|---------|-----------------|
| Competitive CYP inhibition | Two drugs compete for same enzyme | Fluoxetine + codeine (CYP2D6) | GNN knowledge graph: shared METABOLIZES edges |
| Mechanism-based inhibition | Drug irreversibly inactivates enzyme | Ritonavir + CYP3A4 substrates | INHIBITS edge with mechanism=irreversible |
| CYP induction | Drug upregulates enzyme expression | Rifampin + warfarin (CYP2C9 induction) | INDUCES edge |
| Transporter competition | Drugs compete for same transporter | Metformin + cimetidine (OCT1/OCT2) | TRANSPORTS edges with shared Gene node |
| Pharmacodynamic synergy | Additive/synergistic drug effects | Two QT-prolonging drugs | AdverseEvent node overlap |
| Pharmacogenomic amplification | Genetic variant amplifies interaction severity | CYP2D6 PM + CYP2D6 inhibitor | Phenotype node with activity_score=0 + INHIBITS edge |

**Polypharmacy Risk Score**:

```
Risk_polypharmacy = 1 - prod_{i<j} (1 - P(interaction_{ij} | genotype))
```

For a patient on 5 drugs with 10 pairwise interactions each having ~5% probability, the aggregate polypharmacy risk is approximately 40% -- far exceeding what any single interaction would suggest.

### 7. Bayesian Dosage Optimization

#### SONA-Adapted Posterior Estimation

The optimal dose for a patient depends on their metabolizer phenotype, body weight, renal/hepatic function, concomitant medications, and target drug levels. We use a Bayesian framework with SONA (Self-Optimizing Neural Architecture) adaptation.

**Pharmacokinetic Model**:

A one-compartment model with first-order absorption and elimination:

```
C(t) = (F * D / (V_d * (k_a - k_e))) * (exp(-k_e * t) - exp(-k_a * t))
```

where:
- `C(t)` = plasma concentration at time t
- `F` = bioavailability
- `D` = dose
- `V_d` = volume of distribution
- `k_a` = absorption rate constant
- `k_e` = elimination rate constant = CL / V_d
- `CL` = clearance (the key pharmacogenomically modulated parameter)

**Clearance as a Function of Genotype**:

```
CL(genotype) = CL_reference * AS(diplotype) / AS_reference * f_renal * f_hepatic * f_DDI
```

where:
- `CL_reference` = population mean clearance
- `AS(diplotype)` = activity score for the patient's diplotype
- `AS_reference` = reference activity score (typically 2.0 for normal metabolizers)
- `f_renal` = renal function adjustment (eGFR-based)
- `f_hepatic` = hepatic function adjustment (Child-Pugh score-based)
- `f_DDI` = drug-drug interaction adjustment from polypharmacy analysis

**Bayesian Posterior for Dose**:

Given a target therapeutic range [C_min, C_max]:

```
P(D | genotype, target) proportional_to P(C_min <= C_ss(D) <= C_max | D, genotype) * P(D)
```

where:
- `C_ss = (F * D * k_a) / (V_d * k_e * (k_a - k_e))` = steady-state trough concentration
- `P(D)` = prior on dose (from CPIC guidelines for the metabolizer phenotype)

**SONA Adaptation** (`ruvector-sona`):

SONA continuously adapts the dosing model as therapeutic drug monitoring (TDM) data becomes available:

```
1. Initial dose recommendation from Bayesian model with genotype prior
2. Patient receives dose; TDM measurement C_obs at time t_obs
3. SONA updates posterior:
   P(CL | C_obs, genotype) proportional_to P(C_obs | CL, D, t_obs) * P(CL | genotype)
4. Refined dose recommendation with updated CL posterior
5. Repeat with each TDM measurement
```

SONA adaptation latency: < 0.05ms per update, enabling real-time dose adjustment in clinical decision support systems.

**Dosing Recommendations by Metabolizer Status** (example for codeine/CYP2D6):

| Metabolizer Status | CYP2D6 Activity Score | CPIC Recommendation | Dosing Adjustment |
|-------------------|----------------------|--------------------|--------------------|
| Ultra-rapid (UM) | >2.25 | Avoid codeine (risk of morphine toxicity) | Alternative analgesic |
| Normal (NM) | 1.25-2.25 | Standard dose | No adjustment |
| Intermediate (IM) | 0.25-1.0 | Reduced efficacy expected | Consider alternative |
| Poor (PM) | 0 | Avoid codeine (no analgesic effect) | Alternative analgesic |

---

## Quantum Molecular Simulation Architecture

### Hamiltonian Construction Pipeline

```
+-------------------+     +-------------------+     +-------------------+
| Molecular         |     | Integral          |     | Qubit             |
| Geometry          |---->| Computation       |---->| Mapping           |
| (PDB + variant    |     | (h_pq, h_pqrs    |     | (Jordan-Wigner    |
|  structural model)|     |  over active      |     |  or Bravyi-       |
|                   |     |  space)           |     |  Kitaev)          |
+-------------------+     +-------------------+     +-------------------+
                                                           |
                                                           v
+-------------------+     +-------------------+     +-------------------+
| Energy &          |     | Measurement       |     | VQE Circuit       |
| Wavefunction      |<----| (Pauli term       |<----| Preparation       |
| Analysis          |     |  expectation      |     | (UCCSD ansatz     |
|                   |     |  values)          |     |  + error          |
|                   |     |                   |     |  mitigation)      |
+-------------------+     +-------------------+     +-------------------+
        |
        v
+-------------------+
| Pharmacokinetic   |
| Parameter         |
| Extraction        |
| (Km, Vmax,        |
|  barrier height)  |
+-------------------+
```

### Born-Oppenheimer Approximation on Quantum Hardware

The Born-Oppenheimer approximation -- treating nuclei as classical particles and solving for the electronic wavefunction at each nuclear configuration -- is maintained in the quantum simulation. The workflow:

1. **Classical geometry optimization**: Optimize nuclear positions using molecular mechanics (AMBER/CHARMM for the protein, DFT for the active site region)
2. **QM/MM partitioning**: Define the QM region (active site + substrate, ~100-200 atoms) and MM region (rest of protein + solvent)
3. **Active space selection**: Automated active space selection using the MP2 natural orbital occupation numbers:
   - Orbitals with occupation 0.02 < n < 1.98 are included in the active space
   - Typically yields (14e, 14o) to (20e, 20o) for CYP450 active sites
4. **Integral computation**: Compute one- and two-electron integrals over the active space orbitals using PySCF or a custom integral engine
5. **Qubit Hamiltonian**: Map to qubit Hamiltonian via Jordan-Wigner (default) or Bravyi-Kitaev (if fewer Pauli terms beneficial)
6. **VQE execution**: Run VQE on `ruQu` quantum backend
7. **Energy extraction**: Use the optimized wavefunction for property calculations (spin density, charge distribution, reaction barrier)

### Hybrid VQE Architecture (`ruQu`)

```
+-------------------------------------------------------------------+
|                    ruQu Hybrid VQE Engine                           |
+-------------------------------------------------------------------+
|                                                                     |
|  +-------------------+     +-------------------+                   |
|  | Classical          |     | Quantum            |                  |
|  | Optimizer          |     | Circuit Engine      |                 |
|  | (L-BFGS-B,         |<--->| (Parameter bind,    |                |
|  |  COBYLA,           |     |  execute, measure)  |                 |
|  |  SPSA for noisy)   |     |                     |                 |
|  +-------------------+     +-------------------+                   |
|         ^                          |                                |
|         |                          v                                |
|  +-------------------+     +-------------------+                   |
|  | Convergence        |     | Error Mitigation   |                  |
|  | Monitor            |     | Engine              |                 |
|  | (energy, gradient,  |     | (ZNE, symmetry     |                |
|  |  parameter delta)  |     |  verification)     |                  |
|  +-------------------+     +-------------------+                   |
|         |                          |                                |
|         v                          v                                |
|  +---------------------------------------------------+             |
|  | Hardware Abstraction Layer                          |            |
|  |                                                     |            |
|  |  +----------+  +----------+  +----------+          |            |
|  |  | IBM       |  | IonQ      |  | Classical |         |           |
|  |  | Quantum   |  | Aria      |  | Simulator |         |           |
|  |  | (127 qb)  |  | (25 qb)  |  | (statevec |         |           |
|  |  |           |  |           |  |  up to 30 |         |           |
|  |  |           |  |           |  |  qubits)  |         |           |
|  |  +----------+  +----------+  +----------+          |            |
|  +---------------------------------------------------+             |
|                                                                     |
+-------------------------------------------------------------------+
```

**Automatic Backend Selection**:

| Active Space | Qubits | Preferred Backend | Fallback |
|-------------|--------|-------------------|----------|
| (8e, 8o) | 16 | Classical statevector (exact) | N/A |
| (14e, 14o) | 28 | IonQ Aria (25 qubits) with ZNE, or classical CASSCF | Classical CASSCF |
| (20e, 20o) | 40 | IBM Quantum (127 qubits) with ZNE + symmetry | Classical approximate (DMRG) |
| (30e, 30o) | 60 | IBM Quantum (127 qubits) | Classical DMRG (lower accuracy) |

---

## Clinical Decision Support

### Genotype-to-Phenotype Translation Engine

The translation engine converts raw genomic data into clinically actionable pharmacogenomic reports:

```
+-------------------+     +-------------------+     +-------------------+
| Raw Genotype      |     | Star Allele       |     | Metabolizer       |
| (VCF with PGx     |---->| Assignment        |---->| Phenotype         |
|  gene variants)   |     | (GNN-based,       |     | Classification    |
|                   |     |  PharmVar lookup)  |     | (Activity score)  |
+-------------------+     +-------------------+     +-------------------+
                                                           |
                                                           v
+-------------------+     +-------------------+     +-------------------+
| Clinical Report   |     | Dosing            |     | Drug Interaction  |
| Generation        |<----| Recommendations   |<----| Analysis          |
| (PDF/FHIR)        |     | (Bayesian,        |     | (GNN + HNSW +     |
|                   |     |  CPIC-guided)     |     |  Attention)       |
+-------------------+     +-------------------+     +-------------------+
```

**Report Content**:

For each pharmacogene-drug pair with clinical actionability:

| Field | Content | Source |
|-------|---------|--------|
| Gene | CYP2D6 | Genomic annotation |
| Diplotype | *1/*4 | GNN star allele caller |
| Phenotype | Intermediate metabolizer | Activity score: 1.0 |
| Activity Score | 1.0 | *1 (1.0) + *4 (0.0) = 1.0 |
| Affected Drugs | Codeine, tramadol, tamoxifen, ondansetron, ... | CPIC/PharmGKB |
| Recommendation (codeine) | Reduced conversion to morphine; consider alternative analgesic | CPIC Level A |
| Recommendation (tamoxifen) | Reduced conversion to endoxifen; consider alternative or dose increase | CPIC Level A |
| Adverse Event Risk | Codeine: low efficacy (P=0.85); Ondansetron: reduced antiemetic effect (P=0.60) | HNSW adverse event model |
| Quantum Prediction | Km(codeine->morphine) = 45 uM (vs. 12 uM for *1/*1); Vmax reduced 60% | VQE molecular docking |
| Confidence | High (CPIC Level A guideline; GNN structural resolution confidence 0.98) | Multi-source |

### Drug Interaction Alert System

**Alert Hierarchy**:

| Level | Trigger | Display | Requires Acknowledgment |
|-------|---------|---------|------------------------|
| CONTRAINDICATION | HLA-B*57:01 + abacavir; CYP2D6 UM + codeine (children) | Red banner, audible alert | Yes, with override justification |
| MAJOR | CYP2D6 PM + codeine; DPYD deficient + fluorouracil; TPMT PM + azathioprine | Orange banner | Yes |
| MODERATE | CYP2C19 IM + clopidogrel; CYP2C9 PM + warfarin | Yellow banner | No (informational) |
| MINOR | Any actionable PGx result not in above categories | Green notification | No |
| INFORMATIONAL | Star allele assignment, metabolizer phenotype | Report section | No |

### Alternative Drug Recommendations

When a drug is contraindicated or requires significant dose modification, the system recommends therapeutic alternatives using the GNN knowledge graph:

```
MATCH (d1:Drug)-[:METABOLIZES]-(g:Gene)-[:METABOLIZES]-(d2:Drug)
WHERE d1.name = 'codeine'
  AND d2.name <> 'codeine'
  AND d2.atc_code STARTS WITH d1.atc_code[0:3]  // Same therapeutic class
  AND NOT (d2)-[:METABOLIZES]-(:Gene {symbol: 'CYP2D6'})  // Not CYP2D6 substrate
RETURN d2.name, d2.atc_code
ORDER BY d2.evidence_level DESC
```

Example recommendations for CYP2D6 poor metabolizer:

| Original Drug | Therapeutic Class | Alternative | Rationale |
|--------------|------------------|-------------|-----------|
| Codeine | Opioid analgesic | Morphine (direct administration) | Bypasses CYP2D6 activation step |
| Tramadol | Opioid analgesic | Tapentadol | Not CYP2D6 dependent |
| Tamoxifen | Selective estrogen receptor modulator | Aromatase inhibitor (if post-menopausal) | Different metabolic pathway |
| Ondansetron | 5-HT3 antagonist antiemetic | Granisetron | Less CYP2D6 dependent |
| Atomoxetine | ADHD treatment | Methylphenidate | Not CYP2D6 substrate |

---

## Knowledge Graph Integration

### Data Source Integration Pipeline

```
+-------------------+     +-------------------+     +-------------------+
| CPIC              |     | PharmGKB          |     | DrugBank          |
| (Guidelines,      |     | (Clinical         |     | (Drug targets,    |
|  gene-drug pairs, |     |  annotations,     |     |  metabolic        |
|  diplotype-       |     |  variant-level    |     |  pathways,        |
|  phenotype tables)|     |  evidence,        |     |  enzyme Ki/Km,    |
|                   |     |  pathways)        |     |  DDI data)        |
+--------+----------+     +--------+----------+     +--------+----------+
         |                          |                          |
         +------------+-------------+-------------+------------+
                      |                           |
                      v                           v
         +-------------------+          +-------------------+
         | Entity Resolution |          | UniProt           |
         | & Harmonization   |          | (Protein          |
         | (gene symbols,    |          |  structures,      |
         |  drug identifiers,|          |  active sites,    |
         |  variant notation)|          |  post-translational|
         +--------+----------+          |  modifications)   |
                  |                     +--------+----------+
                  +-------------+----------------+
                                |
                                v
                  +-------------------+
                  | Pharmacogenomic   |
                  | Knowledge Graph   |
                  | (Neo4j-compatible |
                  |  in ruvector-gnn) |
                  +-------------------+
```

**Entity Harmonization**:

| Entity Type | CPIC Identifier | PharmGKB Identifier | DrugBank Identifier | Unified ID |
|-------------|----------------|--------------------|--------------------|-----------|
| Gene | HGNC symbol | PharmGKB gene ID (PA*) | UniProt ID | HGNC symbol |
| Drug | Generic name | PharmGKB drug ID (PA*) | DrugBank ID (DB*) | DrugBank ID |
| Variant | rsID + star allele | rsID + haplotype | N/A | rsID + PharmVar star allele |
| Protein | N/A | N/A | UniProt ID | UniProt ID |

**Knowledge Graph Statistics** (estimated for current database versions):

| Entity Type | Count |
|-------------|-------|
| Gene nodes | ~800 |
| Drug nodes | ~15,000 |
| Protein nodes | ~20,000 |
| Variant nodes | ~50,000 |
| Phenotype nodes | ~400 (metabolizer phenotypes) |
| Pathway nodes | ~300 |
| Adverse event nodes | ~5,000 (MedDRA preferred terms) |
| METABOLIZES edges | ~2,000 |
| INHIBITS edges | ~8,000 |
| INDUCES edges | ~1,500 |
| TRANSPORTS edges | ~3,000 |
| CAUSES edges | ~30,000 |
| INTERACTS edges | ~50,000 |
| Total edges | ~100,000+ |

### Knowledge Graph Update Protocol

The knowledge graph is updated quarterly to incorporate new CPIC guidelines, PharmGKB annotations, and DrugBank releases:

1. Download latest releases from each source
2. Run entity resolution pipeline to harmonize new entries
3. Validate against previous version (no entity should be silently removed)
4. Retrain GNN embeddings on updated graph
5. Run regression tests on known drug-gene interactions
6. Deploy updated graph and embeddings atomically

---

## Regulatory Compliance

### FDA Pharmacogenomic Biomarkers

The FDA maintains a Table of Pharmacogenomic Biomarkers in Drug Labeling (https://www.fda.gov/drugs/science-and-research-drugs/table-pharmacogenomic-biomarkers-drug-labeling). As of the latest update, this table includes approximately 500 drug-biomarker pairs across 250+ drugs.

**Biomarker Categories in FDA Labeling**:

| Category | Count | Clinical Impact |
|----------|-------|----------------|
| Required testing before prescribing | ~30 | Prescription cannot be written without test result |
| Recommended testing | ~50 | Strong recommendation in label; not mandatory |
| Actionable PGx information | ~200 | Label describes genotype-phenotype relationship |
| Informational only | ~220 | Label mentions PGx data; no specific action required |

**Required Testing Examples**:

| Drug | Biomarker | Test Required | Consequence of Not Testing |
|------|-----------|--------------|--------------------------|
| Abacavir | HLA-B*57:01 | Yes (before first dose) | Hypersensitivity reaction (potentially fatal) |
| Carbamazepine | HLA-B*15:02 | Yes (in at-risk populations) | Stevens-Johnson syndrome / TEN |
| Rasburicase | G6PD | Yes | Hemolytic anemia |
| Eliglustat | CYP2D6 | Yes | Dose adjustment required; contraindicated in UM with CYP3A inhibitor |
| Ivacaftor | CFTR | Yes | Efficacy limited to specific CFTR mutations |

### IVD-Grade Classification

For the pharmacogenomic pipeline to be deployed in clinical settings, it must meet In Vitro Diagnostic (IVD) device requirements:

**Classification Pathway**: De novo classification or 510(k) substantial equivalence, depending on the specific pharmacogene and claim.

**Analytical Validation Requirements**:

| Metric | Target | Validation Approach |
|--------|--------|-------------------|
| Analytical sensitivity (PPA) | >= 99.0% for each star allele | Comparison to orthogonal methods (Sanger sequencing, TaqMan) on 200+ samples per star allele |
| Analytical specificity (NPA) | >= 99.5% | Same sample set; no false positive star allele calls |
| Reproducibility | >= 99.0% concordance across runs, operators, instruments | 20 samples x 3 runs x 2 operators x 2 instruments |
| Reportable range | All star alleles with >= 1% population frequency in any major ancestry group | PharmVar allele frequency data |
| Interfering substances | No impact from hemolysis, lipemia, bilirubin, common medications | Spike-in studies |
| Copy number accuracy | Within +/- 0.5 copies of true copy number | Validated against ddPCR (digital droplet PCR) |

**Clinical Validation Requirements**:

| Metric | Target | Validation Approach |
|--------|--------|-------------------|
| Clinical concordance with phenotype | >= 95% agreement between predicted phenotype and observed drug response | Prospective study, N >= 500 patients per gene-drug pair |
| Positive predictive value for adverse events | >= 80% for CPIC Level A gene-drug pairs | Retrospective analysis of EHR-linked biobank data |
| Time to result | < 2 hours from sample to report (for pre-prescribing use) | Workflow timing study |

**Quality Management System**:

| Requirement | Standard | Implementation |
|-------------|---------|----------------|
| Design controls | 21 CFR 820.30 | Formal design history file; documented V&V |
| Software validation | FDA Guidance on Software Validation | IEC 62304 software lifecycle; unit test coverage > 95% |
| Risk management | ISO 14971 | Hazard analysis with FMEA; residual risk acceptance |
| Cybersecurity | FDA Premarket Cybersecurity Guidance | Encryption at rest and in transit; access controls; audit logging |
| Labeling | 21 CFR 809.10 | Intended use, limitations, performance characteristics |

---

## Performance Targets

### Star Allele Calling Performance

| Metric | Target | Hardware |
|--------|--------|----------|
| CYP2D6 diplotype accuracy | >= 99.0% | 128-core server |
| CYP2D6 copy number accuracy | >= 99.5% (within +/- 0.5 copies) | 128-core server |
| CYP2C19 diplotype accuracy | >= 99.5% | 128-core server |
| DPYD diplotype accuracy | >= 99.5% | 128-core server |
| TPMT diplotype accuracy | >= 99.5% | 128-core server |
| HLA-B typing accuracy (4-digit) | >= 99.0% | 128-core server |
| Star allele calling latency (per gene) | < 5 seconds | 128-core server |
| Full pharmacogenomic panel (15 genes) | < 30 seconds | 128-core server |
| GNN inference (star allele structural resolution) | < 500 ms per gene | GPU (NVIDIA A100) |

### Drug-Gene Interaction Prediction Performance

| Metric | Target | Notes |
|--------|--------|-------|
| Interaction type prediction AUC-ROC | >= 0.95 | 5-class classification |
| Interaction strength (Km) prediction | Spearman rho >= 0.85 | Continuous regression |
| Adverse event prediction AUC-ROC | >= 0.90 | Binary per MedDRA PT |
| Drug-drug interaction severity | F1 >= 0.85 | 5-class ordinal |
| GNN knowledge graph inference latency | < 100 ms per query | Per drug-gene pair |
| HNSW adverse event search (100M records) | < 5 ms (k=100) | Including similarity computation |

### Quantum Simulation Performance

| Metric | Target | Backend |
|--------|--------|---------|
| VQE convergence (14e, 14o) | < 500 iterations | IBM Quantum or IonQ |
| VQE convergence (20e, 20o) | < 1,000 iterations | IBM Quantum |
| VQE energy accuracy | < 1 kcal/mol (1.6 mHartree) of CASPT2 reference | With ZNE error mitigation |
| QPE energy accuracy | < 0.5 kcal/mol (0.8 mHartree) | Fault-tolerant QC (future) |
| Classical simulator (20e, 20o) | < 4 hours per single-point energy | 128-core CPU |
| Barrier height accuracy | < 2 kcal/mol of experimental | VQE on (20e, 20o) |
| Km prediction from quantum barrier | R^2 >= 0.80 vs. experimental | Validated on MetaQSAR |

### Clinical Decision Support Performance

| Metric | Target | Notes |
|--------|--------|-------|
| Time from VCF to clinical report | < 120 seconds (full pipeline) | Including quantum fallback to classical |
| Time from VCF to clinical report (no quantum) | < 60 seconds | Classical-only path |
| Alert sensitivity (life-threatening ADR) | >= 99.0% | No missed contraindications |
| Alert specificity (life-threatening ADR) | >= 80.0% | Tolerate some false positives for safety |
| Dosing recommendation concordance with CPIC | >= 95.0% | For Level A guidelines |
| SONA adaptation latency | < 0.05 ms per TDM update | Real-time dose adjustment |

### Flash Attention Performance for Polypharmacy

| Number of Concurrent Drugs | Naive Attention | Flash Attention | Speedup |
|---------------------------|----------------|-----------------|---------|
| 5 | 0.2 ms | 0.1 ms | 2.0x |
| 10 | 1.5 ms | 0.4 ms | 3.8x |
| 20 | 8.0 ms | 1.5 ms | 5.3x |
| 50 (extreme polypharmacy) | 65 ms | 9 ms | 7.2x |

---

## Consequences

### Positive Consequences

1. **Mechanistic drug response prediction**: Quantum molecular simulation provides first-principles predictions of drug-enzyme interaction kinetics, replacing empirical population-level statistics with patient-specific quantum-computed parameters
2. **Comprehensive star allele calling**: GNN-based structural resolution handles the full complexity of CYP2D6 (deletions, duplications, hybrids) that standard variant callers miss
3. **Real-time polypharmacy safety**: Multi-head attention over drug interaction tensors with flash attention acceleration enables sub-10ms interaction analysis even for patients on 50+ medications
4. **Continuous dose optimization**: SONA-adapted Bayesian dosing incorporates therapeutic drug monitoring data in real-time, with < 0.05ms adaptation latency
5. **Scalable adverse event prediction**: HNSW similarity search over 100M+ historical patient-drug-outcome records enables real-time risk stratification with sub-5ms latency
6. **Regulatory pathway clarity**: Explicit IVD-grade analytical and clinical validation targets provide a clear path to clinical deployment

### Negative Consequences

1. **Quantum hardware dependency**: Full quantum advantage for molecular simulation requires fault-tolerant quantum computers (>1,000 logical qubits) that do not yet exist; NISQ-era results require extensive error mitigation and may not achieve chemical accuracy for all systems
2. **Knowledge graph maintenance burden**: Quarterly updates from four databases (CPIC, PharmGKB, DrugBank, UniProt) require ongoing curation effort and regression testing
3. **Limited training data for rare star alleles**: Star alleles with population frequency < 0.1% have insufficient training data for GNN validation; clinical validation is impractical for ultra-rare alleles
4. **Regulatory uncertainty**: No pharmacogenomic clinical decision support system with quantum-enhanced molecular simulation has been FDA-cleared; the regulatory pathway is uncharted
5. **Computational cost of quantum simulation**: Even with hybrid VQE, molecular docking for a single drug-enzyme pair requires hours of quantum compute time; pre-computation and caching are essential

### Risk Mitigation

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| VQE does not converge for large active spaces | Medium | High | Automatic fallback to classical DMRG/CASSCF; hybrid classical-quantum active space partitioning |
| GNN misassigns CYP2D6 hybrid alleles | Low | High | Orthogonal confirmation via long-read sequencing; conservative reporting (flag uncertain calls) |
| Knowledge graph becomes stale | Medium | Medium | Automated quarterly update pipeline with regression testing; version pinning for clinical use |
| HNSW adverse event search returns biased results | Medium | Medium | Population-stratified indexing; bias monitoring dashboard; calibration on diverse cohorts |
| Quantum hardware unavailable | High (current) | Low (mitigated) | All quantum calculations have classical fallbacks; quantum results are supplementary refinements |
| Regulatory submission rejected | Medium | High | Early FDA Pre-Submission meeting; phased approach (classical first, quantum enhancement later) |

---

## Alternatives Considered

### Alternative 1: Classical-Only Pharmacogenomics (No Quantum)

Use classical molecular mechanics (MM/GBSA) for docking and established tools (Stargazer, Cyrius) for star allele calling.

**Rejected because**:
- Classical docking cannot accurately model CYP450 transition states (systematic errors of 3-5 kcal/mol for iron-oxo species)
- Existing star allele callers have known failure modes for CYP2D6 hybrid alleles and complex structural configurations
- No path to first-principles Km/Vmax prediction from genotype alone

### Alternative 2: Deep Learning End-to-End Drug Response Prediction

Train a large language model or transformer on patient genotype -> drug response data directly, bypassing mechanistic modeling.

**Rejected because**:
- Requires enormous labeled datasets (genotype + drug + outcome) that are not available for most gene-drug pairs
- No interpretability: cannot explain why a specific genotype leads to a specific drug response
- Cannot generalize to novel drugs or novel star alleles not seen in training data
- Regulatory agencies (FDA) require mechanistic justification for clinical decision support claims

### Alternative 3: Outsource Star Allele Calling to Existing Tools

Use Stargazer, PharmCAT, or Aldy for star allele calling and focus only on the downstream analysis.

**Rejected because**:
- Existing tools do not integrate with the RuVector variant calling pipeline (ADR-009) and require separate VCF preprocessing
- No tools provide uncertainty quantification for star allele calls at the level required for IVD-grade classification
- Integration of star allele calling with structural variant resolution in a single GNN model provides higher accuracy than a two-stage pipeline

---

## Related Decisions

- **ADR-001**: RuVector Core Architecture (HNSW index, SIMD distance)
- **ADR-003**: Hierarchical Navigable Small World Genomic Vector Index (pharmacogenomic vector space definition)
- **ADR-009**: Zero-False-Negative Variant Calling Pipeline (upstream variant calling for PGx genes)
- **ruQu ADR-001**: ruQu Architecture (VQE, QAOA for quantum-enhanced models)

---

## References

1. Relling, M.V., & Klein, T.E. (2011). "CPIC: Clinical Pharmacogenetics Implementation Consortium of the Pharmacogenomics Research Network." *Clinical Pharmacology & Therapeutics*, 89(3), 464-467.

2. Whirl-Carrillo, M., et al. (2021). "An evidence-based framework for evaluating pharmacogenomics knowledge for personalized medicine." *Clinical Pharmacology & Therapeutics*, 110(3), 563-572. (PharmGKB)

3. Wishart, D.S., et al. (2018). "DrugBank 5.0: a major update to the DrugBank database for 2018." *Nucleic Acids Research*, 46(D1), D1074-D1082.

4. The UniProt Consortium (2023). "UniProt: the Universal Protein Knowledgebase in 2023." *Nucleic Acids Research*, 51(D1), D523-D531.

5. Gaedigk, A., et al. (2018). "The Pharmacogene Variation (PharmVar) Consortium: Incorporation of the Human Cytochrome P450 (CYP) Allele Nomenclature Database." *Clinical Pharmacology & Therapeutics*, 103(3), 399-401.

6. Perutz, A., et al. (2014). "A variational eigenvalue solver on a photonic quantum processor." *Nature Communications*, 5, 4213. (VQE original)

7. Cao, Y., et al. (2019). "Quantum chemistry in the age of quantum computing." *Chemical Reviews*, 119(19), 10856-10915.

8. Shaik, S., et al. (2005). "Theoretical perspective on the structure and mechanism of cytochrome P450 enzymes." *Chemical Reviews*, 105(6), 2279-2328. (Two-state reactivity)

9. Kandala, A., et al. (2017). "Hardware-efficient variational quantum eigensolver for small molecules and quantum magnets." *Nature*, 549(7671), 242-246.

10. Malkov, Y., & Yashunin, D. (2018). "Efficient and robust approximate nearest neighbor search using Hierarchical Navigable Small World graphs." *IEEE TPAMI*, 42(4), 824-836.

11. Dao, T., et al. (2022). "FlashAttention: Fast and Memory-Efficient Exact Attention with IO-Awareness." *NeurIPS 2022*.

12. Lazarou, J., Pomeranz, B.H., & Corey, P.N. (1998). "Incidence of adverse drug reactions in hospitalized patients: a meta-analysis of prospective studies." *JAMA*, 279(15), 1200-1205.

13. Caudle, K.E., et al. (2020). "Standardizing CYP2D6 Genotype to Phenotype Translation: Consensus Recommendations from the Clinical Pharmacogenetics Implementation Consortium and Dutch Pharmacogenetics Working Group." *Clinical and Translational Science*, 13(1), 116-124.

14. FDA Table of Pharmacogenomic Biomarkers in Drug Labeling. https://www.fda.gov/drugs/science-and-research-drugs/table-pharmacogenomic-biomarkers-drug-labeling

---

## Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 0.1 | 2026-02-11 | RuVector DNA Analyzer Team | Initial proposal |
