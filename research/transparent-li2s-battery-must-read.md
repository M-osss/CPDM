# Must-read documents: transparent lithium sulfide battery research

Compiled 2026-09-11.

There is no published, fully transparent Li2S cell. The field is the intersection of three literatures: (1) transparent / semitransparent thin-film batteries, (2) Li2S as a prelithiated sulfur cathode, (3) optical and electronic structure of antifluorite Li2S. Pure crystalline Li2S is a wide-gap insulator (calculated indirect gap ~3.5–3.7 eV) and is therefore optically transparent in the visible if chemically pure. Off-white or yellow powder is a polysulfide impurity signature, not an intrinsic optical property.

Theoretical Li2S capacity is 1166 mAh g−1. Density 1.66 g cm−3. Melting point 938 °C. Electronic conductivity is negligible (~10−24 S m−1 in the ideal crystal). First-charge activation typically requires ~3.5 V vs Li/Li+. Moisture converts Li2S to LiOH + H2S.

A transparent Li2S device must therefore solve conductivity without destroying visible transmittance: grid patterning of opaque conductors, ultrathin vapor-deposited Li2S on ITO/FTO/graphene, or a solid-state stack that keeps carbon/metal below the eye-resolution limit.

---

## Read first (in this order)

| # | Document | Why it is first |
|---|----------|-----------------|
| 1 | Yang et al., *PNAS* 2011 | Defines the transparency–capacity dilemma and the sub-eye-resolution grid architecture. 60% transmittance at 10 Wh L−1. |
| 2 | Oukassi et al., *ACS AMI* 2019 | First all-inorganic transparent thin-film Li-ion cell (LiCoO2/LiPON/Si). 25–60% T. Industrial PVD + lithography path. |
| 3 | Yang et al., *JACS* 2012 | Li2S first-charge overpotential (~1 V extra). Activation by raising the charge cutoff; thereafter the barrier vanishes. |
| 4 | Ting et al., *ACS Omega* 2022 | Single Li2S-cathode review covering synthesis, hosts, additives, and the activation-barrier literature. |
| 5 | Klein & Manthiram, *JACS* 2017 | Sputtered Li2S films: as-deposited polymer-like sulfide is electrochemically active; annealed crystalline Li2S is inactive. Direct thin-film process document. |
| 6 | Deng et al., *Nano-Micro Lett.* 2023 | First demonstrated all-solid-state thin-film Li–S cell: VGs-Li2S / LiPON / Li. 81% retention over 3000 cycles. Closest existing hardware to a transparent Li2S stack. |
| 7 | Bates et al., *J. Power Sources* 1993 | Origin of LiPON. ~2×10−6 S cm−1, Li-metal stable, optically usable as a thin amorphous film. |
| 8 | Eithiraj et al., *phys. status solidi (b)* 2007 | Electronic structure of antifluorite Li2S. Indirect gap. Ground-state lattice and bulk modulus. |
| 9 | Wehner, Mittal, Liu, Niederberger, *ACS Cent. Sci.* 2021 | Transparent-battery materials review: current collectors, polymer vs inorganic electrolytes, transmittance/resistance trade-offs. |
| 10 | Kim, Su, Zhong et al., *Nat. Chem. Eng.* 2024 | Solid-state sulfur redox as a transport–kinetics problem (Damköhler number). Needed once the cell is all-solid-state. |

---

## A. Transparent battery architecture

These papers establish how to make a battery look transparent. None of them use Li2S. They supply the optical design rules.

1. **Yang, Y.; Jeong, S.; Hu, L.; Wu, H.; Lee, S. W.; Cui, Y.** Transparent lithium-ion batteries. *Proc. Natl. Acad. Sci. U.S.A.* **2011**, *108*, 13013–13018.  
   https://doi.org/10.1073/pnas.1102873108  
   Open: https://pmc.ncbi.nlm.nih.gov/articles/PMC3156205/  
   Grid electrodes with feature size below human-eye resolution. Aligning stacked grids raises energy without a linear loss of transmittance. Flexible PDMS substrate. Also usable as an in-situ Raman window.

2. **Oukassi, S.; Baggetto, L.; Dubarry, C.; Le Van-Jodin, L.; Poncet, S.; Salot, R.** Transparent thin film solid-state lithium ion batteries. *ACS Appl. Mater. Interfaces* **2019**, *11*, 683–690.  
   https://doi.org/10.1021/acsami.8b16364  
   LiCoO2/LiPON/Si on glass. Photolithography + etch to sub-eye grids. UV–vis T up to 60%. Discharge ~0.15 mAh at C/2 (3–4.2 V) for the highest-T variants. Architecture-parameter sweep of T vs capacity. This is the fabrication template for an inorganic transparent Li2S cell: replace LCO with sputtered/ALD Li2S.

3. **Oukassi, S. et al.** Transparent thin film lithium ion batteries. *ECS Meeting Abstracts* **2017**, MA2017-02, 56.  
   https://doi.org/10.1149/MA2017-02/1/56  
   8-inch glass process note. LCO up to 10 µm, Si ~0.1 µm. T range 25–60%. Read after the 2019 AMI paper.

4. **Pat, S. et al.** The microstructural, surface, optical and electrochemical impedance spectroscopic study of the semitransparent all-solid-state thin film battery. *Mater. Res. Express* **2019**, *6*, 015503.  
   https://doi.org/10.1088/2053-1591/aae4aa  
   Semitransparent polycrystalline TFB. ~8.3 µAh cm−2, ~100% CE, 100 cycles. Optical reflectance of the full stack vs bare glass. Photovoltaic-window application framing.

5. **Wehner, L. A.; Mittal, N.; Liu, T.; Niederberger, M.** Multifunctional batteries: flexible, transient, and transparent. *ACS Cent. Sci.* **2021**, *7*, 231–244.  
   https://doi.org/10.1021/acscentsci.0c01318  
   Open: https://pmc.ncbi.nlm.nih.gov/articles/PMC7908028/  
   Compact survey of transparent current collectors: ITO/FTO, metal nanowires (~24.5 Ω sq−1 at 71% T), CNT (57 Ω sq−1 at 90% T), graphene (350 Ω sq−1 at 90% T), PEDOT:PSS (260 Ω sq−1 at 95% T). Transparent electrolytes: Al2O3-doped LLZO (~30% T, 9.9×10−4 S cm−1), PVDF-HFP / PEO / PMMA gels. Use this to pick the collector before depositing Li2S. ETH funding line in the paper is literally “Towards transparent lithium ion batteries.”

6. **Jia, B.; Zhang, C.; Liu, M.; Li, Z.; Wang, J.; Zhong, L.; Han, C.; Qin, M.; Huang, X.** Integration of microbattery with thin-film electronics for constructing an integrated transparent microsystem based on InGaZnO. *Nat. Commun.* **2023**, *14*, 5330.  
   https://doi.org/10.1038/s41467-023-41181-1  
   Open: https://pmc.ncbi.nlm.nih.gov/articles/PMC10474284/  
   Transparent IGZO anode LIB (~9.8 µAh cm−2) monolithically integrated with TFT and photodetector. Shows what a transparent Li2S microbattery would have to mate with.

7. **Lee, S. et al.** Solid state thin electrolyte to overcome transparency-capacity dilemma of transparent supercapacitor. *Sci. Rep.* **2022**, *12*, 15854.  
   https://doi.org/10.1038/s41598-022-19933-8  
   ITO/LiCoO2/LiPON/WO3/ITO. >60% T. LiPON thickness as the optical–ionic lever. Same PVD stack language as a transparent Li2S cell, with WO3 as a transparent intercalation counter-electrode.

Design rules extracted from A:

- Opaque active material is allowed if lateral feature width ≲ 35–50 µm (below eye resolution).
- Stacked, registered grids raise energy at nearly constant T.
- Continuous thin films need every layer (collector, cathode, electrolyte, anode) to be wide-gap or sub-100 nm.
- LiPON is the default transparent inorganic electrolyte for PVD stacks.
- ITO/FTO on glass or PET is the default collector; graphene/AgNW if flexibility is required.

---

## B. Li2S cathode electrochemistry

Li2S is already fully lithiated. It pairs with Li-free anodes (Si, Sn, graphite, anode-less Cu) and starts at maximum volume, so charge-induced expansion is smaller than in S8 cells. The penalties: insulator, moisture-sensitive, high first-charge overpotential, shuttle still present in liquid electrolytes.

### Reviews (read these, then pick primary papers)

8. **Ting, L. K. J.; Gao, Y.; Wang, H.; Wang, T.; Sun, J.; Wang, J.** Lithium sulfide batteries: addressing the kinetic barriers and high first charge overpotential. *ACS Omega* **2022**, *7*, 40682–40700.  
   https://doi.org/10.1021/acsomega.2c05477  
   Open: https://pmc.ncbi.nlm.nih.gov/articles/PMC9670706/  
   Best single Li2S-cathode review. Fundamentals, S vs Li2S comparison, nanoparticle synthesis, conductive matrices, redox mediators, and the Yang vs Zhang activation-mechanism dispute.

9. **Kaiser, M. R.; Han, Z.; Liang, J.; Dou, S. X.; Wang, J.** Lithium sulfide-based cathode for lithium-ion/sulfur battery: recent progress and challenges. *Energy Storage Mater.* **2019**, *19*, 1–15.  
   https://doi.org/10.1016/j.ensm.2019.04.001  
   Full-cell focus: Li-metal-free anodes and electrolytes for Li2S.

10. **Gao, G.; Yang, X.; Bi, J.; Guan, W.; Du, Z.; Ai, W.** Advanced engineering strategies for Li2S cathodes in lithium–sulfur batteries. *J. Mater. Chem. A* **2023**, *11*, 26318–26339.  
    https://doi.org/10.1039/D3TA06057H  
    Polycrystalline Li2S and artificial CEI — topics most reviews skip.

11. **Park, H. et al.** Designing with Li2S in lithium–sulfur batteries: from fundamental chemistry to practical architectures. *Small* **2026**, *22*, e2513644.  
    https://doi.org/10.1002/smll.202513644  
    Open PDF: https://research.chalmers.se/publication/550494/file/550494_Fulltext.pdf  
    2026 perspective: Li2S as a platform for anode-free and solid-state cells. Atomic catalysis, mesoscale architecture, ML-assisted discovery.

12. **Recent progress in host and electrolyte engineering towards Li2S cathode for lithium-sulfur battery.** *Nano Res. Energy* **2025**.  
    https://doi.org/10.26599/NRE.2025.9120194  
    Host + electrolyte only. Use as a 2025 update, not as a first read.

### Li–S system reviews (needed for shuttle, E/S ratio, anode)

13. **Ji, X.; Lee, K. T.; Nazar, L. F.** A highly ordered nanostructured carbon–sulphur cathode for lithium–sulphur batteries. *Nat. Mater.* **2009**, *8*, 500–506.  
    https://doi.org/10.1038/nmat2460  
    CMK-3/S. Origin of the confined-sulfur cathode. 1320 mAh g−1. Every later Li2S–carbon composite is downstream of this paper.

14. **Manthiram, A.; Fu, Y.; Chung, S.-H.; Zu, C.; Su, Y.-S.** Rechargeable lithium–sulfur batteries. *Chem. Rev.* **2014**, *114*, 11751–11787.  
    https://doi.org/10.1021/cr500062v  
    Canonical Li–S review. Electrochemistry, shuttle, anode, electrolyte.

15. **Seh, Z. W.; Sun, Y.; Zhang, Q.; Cui, Y.** Designing high-energy lithium–sulfur batteries. *Chem. Soc. Rev.* **2016**, *45*, 5605–5634.  
    https://doi.org/10.1039/C5CS00410A  
    Materials-design counterpart to Manthiram 2014.

16. **Pang, Q.; Liang, X.; Kwok, C. Y.; Nazar, L. F.** Advances in lithium–sulfur batteries based on multifunctional cathodes and electrolytes. *Nat. Energy* **2016**, *1*, 16132.  
    https://doi.org/10.1038/nenergy2016132  
    Polar hosts and chemical LiPS binding. Sets the language for Li2S deposition morphology.

17. **Ould Ely, T.; Kamzabek, D.; Chakraborty, D.; Doherty, M. F.** Lithium–sulfur batteries: state of the art and future directions. *ACS Appl. Energy Mater.* **2018**, *1*, 1783–1814.  
    https://doi.org/10.1021/acsaem.7b00153  
    Practical metrics: S loading ≥4–6 mg cm−2, S fraction ≥70%, utilization ≥80%, E/S ≤ 3 µL mg−1.

17b. **Zhou, G.; Chen, H.; Cui, Y.** Formulating energy density for designing practical lithium–sulfur batteries. *Nat. Energy* **2022**, *7*, 312–319.  
    https://doi.org/10.1038/s41560-022-01001-0  
    The energy-density accounting paper. Required before claiming any transparent Li2S cell is “high energy”; areal loading and inactive-mass fractions dominate.

### Primary Li2S papers

18. **Yang, Y.; Zheng, G.; Misra, S.; Nelson, J.; Toney, M. F.; Cui, Y.** High-capacity micrometer-sized Li2S particles as cathode materials for advanced rechargeable lithium-ion batteries. *J. Am. Chem. Soc.* **2012**, *134*, 15387–15394.  
    https://doi.org/10.1021/ja3052206  
    Open PDF: http://web.stanford.edu/group/cui_group/papers/Yang_JACS_2012.pdf  
    First-charge kinetic model: insulating Li2S → high charge-transfer resistance → ~3.5 V needed to nucleate soluble polysulfides, which then mediate the rest of the conversion. After activation, overpotential disappears. Micron Li2S, not nano, for tap density.

19. **Zhang, L.; Sun, D.; Feng, J.; Cairns, E. J.; Guo, J.** Revealing the electrochemical charging mechanism of nanosized Li2S by in situ and operando X-ray absorption spectroscopy. *Nano Lett.* **2017**, *17*, 5084–5091.  
    https://doi.org/10.1021/acs.nanolett.7b02381  
    Competing mechanism: first charge is a solid–solid Li2S → S conversion; high overpotential is Li–S bond breaking, not polysulfide nucleation. Read against Yang 2012. Both agree that residual LiPS after the first cycle collapses the barrier.

20. **Zu, C.; Klein, M. J.; Manthiram, A.** Activated Li2S as a high-performance cathode for rechargeable lithium–sulfur batteries. *J. Phys. Chem. Lett.* **2014**, *5*, 3986–3991.  
    https://doi.org/10.1021/jz5021108  
    Chemical/thermal activation of commercial Li2S.

21. **Ye, H.; Li, M.; Liu, T.; Li, Y.; Lu, J.** Activating Li2S as the lithium containing cathode in lithium–sulfur batteries. *ACS Energy Lett.* **2020**, *5*, 2234–2245.  
    https://doi.org/10.1021/acsenergylett.0c00936  
    Concise activation-strategy map (mediators, catalysts, nanostructuring).

22. **Tan, G. et al.** Burning lithium in CS2 for high-performing compact Li2S–graphene nanocapsules for Li–S batteries. *Nat. Energy* **2017**, *2*, 17090.  
    https://doi.org/10.1038/nenergy.2017.90  
    Compact Li2S synthesis. Commentary: Li, Y.; Chen, F. *Nat. Energy* **2017**, *2*, 17096.

23. **Seh, Z. W. et al.** Sulphur–TiO2 yolk–shell nanoarchitecture with internal void space for long-cycle lithium–sulphur batteries. *Nat. Commun.* **2013**, *4*, 1331.  
    https://doi.org/10.1038/ncomms2333  
    Void-space architecture that later transfers to Li2S yolks.

24. **Zhou, W.; Yu, Y.; Chen, H.; DiSalvo, F. J.; Abruña, H. D.** Yolk–shell structure of polyaniline-coated sulfur for lithium–sulfur batteries. *J. Am. Chem. Soc.* **2013**, *135*, 16736–16743.  
    https://doi.org/10.1021/ja409508q  
    Conducting-polymer coating template later used on Li2S.

---

## C. Thin-film Li2S and all-solid-state Li–S

This is the hardware path for a transparent device. Liquid cells with carbon black are optically dead. Vapor-phase Li2S on a TCO, sealed by LiPON, is the realistic stack.

25. **Klein, M. J.; Veith, G. M.; Manthiram, A.** Chemistry of sputter-deposited lithium sulfide films. *J. Am. Chem. Soc.* **2017**, *139*, 9229–9237.  
    https://doi.org/10.1021/jacs.7b03379  
    Custom RF sputter tool. Thickness from a few nm to several µm at >2 nm min−1. As-deposited: polymer-like Li2S chains, full theoretical capacity on first charge, solid-state charge process. Annealed crystalline Li2S: electrochemically inactive. Plasma chemistry produces nonstoichiometry and depth inhomogeneity. This is the process paper for a transparent Li2S cathode.

26. **Klein, M. J.** Understanding the electrochemistry and reaction mechanisms of solid-state sulfides with application to the lithium-sulfur battery system. Ph.D. dissertation, University of Texas at Austin, 2016.  
    https://doi.org/10.15781/t2x63bb46  
    Open: https://repositories.lib.utexas.edu/ (search title).  
    Bandgap 3.5–3.7 eV. Conductivity ~1.9×10−24 S m−1. XPS (Li2S-type vs terminal S), Raman, UV–vis of extracted polysulfides from unannealed vs 600 °C films. Yellow vs white as a purity diagnostic. Required lab-manual companion to the 2017 JACS paper.

27. **Meng, X.; Comstock, D. J.; Fister, T. T.; Elam, J. W.** Vapor-phase atomic-controllable growth of amorphous Li2S for high-performance lithium–sulfur batteries. *ACS Nano* **2014**, *8*, 10963–10972.  
    https://doi.org/10.1021/nn505480w  
    ALD Li2S. ~800 mAh g−1, ~100% CE. Atomic thickness control — the other vapor-phase route besides sputtering. Better for conformal coating of 3D transparent scaffolds.

28. **Deng, R.; Ke, B.; Xie, Y.; Cheng, S.; Zhang, C.; Zhang, H.; Lu, B.; Wang, X.** All-solid-state thin-film lithium-sulfur batteries. *Nano-Micro Lett.* **2023**, *15*, 73.  
    https://doi.org/10.1007/s40820-023-01064-y  
    Open access. VGs-Li2S / LiPON / Li. Li2S chosen over S8 because Tm = 938 °C and it survives sputtering. 81% retention / 3000 cycles with Li-foil anode. Evaporated-Li anode: 500 cycles, 99.71% CE. 20.47 µAh cm−2 at 60 °C. Graphene host is optically absorbing — replace VGs with a grid or a TCO if transparency is the goal. This is the existence proof that a Li2S/LiPON thin-film cell cycles.

29. **Kim, J. T.; Su, H.; Zhong, Y.; Wang, C.; Wu, H.; Zhao, D.; Wang, C.; Sun, X.; Li, Y.** All-solid-state lithium–sulfur batteries through a reaction engineering lens. *Nat. Chem. Eng.* **2024**, *1*, 400–410.  
    https://doi.org/10.1038/s44286-024-00079-5  
    Solid-state S redox ≠ liquid-state. Mass transport, kinetics, thermodynamics. Damköhler number as the design dimensionless group. Cryo-EM as the characterization tool. Read before designing a solid transparent Li2S cathode.

30. **Wang, D.; Jhang, L.-J.; Kou, R. et al.** Realizing high-capacity all-solid-state lithium-sulfur batteries using a low-density inorganic solid-state electrolyte. *Nat. Commun.* **2023**, *14*, 1895.  
    https://doi.org/10.1038/s41467-023-37564-z  
    Liquid-phase Li3PS4–2LiBH4, density 1.491 g cm−3, 6.0 mS cm−1, ~500 nm particles. 60 wt% S, 1144.6 mAh g−1 at 60 °C. Bulk pellet, not transparent, but the low-density SE argument (more SE volume at fixed mass → fewer isolated S domains) transfers to thin-film cathode formulation.

31. **Su, Y. et al.** Progress and prospects of inorganic solid-state electrolyte-based all-solid-state Li–S batteries. *Adv. Sustainable Syst.* **2025**, *9*, 2400555.  
    https://doi.org/10.1002/adsu.202400555  
    Oxide vs sulfide SE, cathode microstructure, Li-anode interlayers, in-situ characterization.

32. **Liang, F. et al.** Insight into all-solid-state Li–S batteries: challenges, advances, and engineering design. *Adv. Energy Mater.* **2024**, *14*, 2401959.  
    https://doi.org/10.1002/aenm.202401959  
    Engineering-design companion to Su 2025.

33. **Nagao, M.; Hayashi, A.; Tatsumisago, M.** High-capacity Li2S–nanocarbon composite electrode for all-solid-state rechargeable lithium batteries. *J. Mater. Chem.* **2012**, *22*, 10015–10020.  
    https://doi.org/10.1039/c2jm16802b  
    Early ASSB Li2S–C composite. Historical baseline.

---

## D. Optical and electronic structure of Li2S

Transparency of a Li2S film is a materials-purity and defect problem, not a dye problem. Wide gap ⇒ visible transparency. S vacancies, polysulfide residues, and carbon hosts destroy it.

34. **Eithiraj, R. D.; Jaiganesh, G.; Kalpana, G.; Rajagopalan, M.** First-principles study of electronic structure and ground-state properties of alkali-metal sulfides — Li2S, Na2S, K2S and Rb2S. *phys. status solidi (b)* **2007**, *244*, 1337–1346.  
    https://doi.org/10.1002/pssb.200642506  
    Antifluorite (anti-CaF2). Li2S is an indirect-gap semiconductor.

35. **Ekuma, C. E.; Jarrell, M.; Moreno, J.; Bagayoko, D.** Ab initio prediction of electronic, transport and bulk properties of Li2S. *Int. J. Mod. Phys. B* **2015**, *29*, 1542006.  
    https://doi.org/10.1142/S0217979215420060  
    Indirect Γ→X gap 3.723 eV at the experimental lattice constant a = 5.689 Å. Effective masses, bulk modulus. LDA with BZW-EF basis, so the gap is not the usual DFT underestimate.

36. **Materials Project mp-1153** (cubic Fm-3m Li2S) and **mp-1125** (Pnma Li2S).  
    https://next-gen.materialsproject.org/materials/mp-1153  
    Computed structures, XRD, band gap (GGA underestimates; treat as a lower bound), elastic tensors. Use as the CIF/POSCAR source, not as the optical ground truth.

37. **Bertheville, B.; Bill, H.; Hagemann, H.** Experimental Raman scattering investigation of phonon anharmonicity effects in Li2S. *J. Phys.: Condens. Matter* **1998**, *10*, 2155–2169.  
    https://doi.org/10.1088/0953-8984/10/9/018  
    Open PDF: https://www.unige.ch/sciences/chifi/publis/refs_pdf/ref00334.pdf  
    Single-crystal antifluorite Li2S. T2g Raman mode 372.6 cm−1 at 295 K, FWHM 10.8 cm−1. Lattice-constant vs T. The fingerprint for “is this film crystalline Li2S?”

38. **See, K. A. et al.** Elucidating the electrochemical activity of electrolyte-insoluble polysulfide species in lithium-sulfur batteries. *J. Electrochem. Soc.* **2016**.  
    https://doi.org/10.1149/2.0051610jes  
    XRD of antifluorite Li2S (JCPDS 23-0369). Explicit statement: a 3.5–3.7 eV gap implies a white ionic solid; yellow/off-white commercial Li2S is a polysulfide impurity marker. Li2S2 calculated gap ~1.8 eV (visible absorption). Use this as the optical-purity criterion for a transparent cathode.

39. **Li, X. et al.** Insight into sulfur vacancy-induced insulator to metal transition of Li2S. *Funct. Mater. Lett.* **2017**, *10*, 1750067.  
    https://doi.org/10.1142/S1793604717500679  
    S vacancy collapses the gap to ~0.7 eV and forms Li–Li metallic bonds. Defects that raise electronic conductivity also kill transparency. This is the central materials trade-off for a transparent Li2S electrode.

40. **First-principles investigation of structural, electronic, optical, and thermoelectric properties of Li2X (X = Te, Se, S).** *Results Phys.* **2025**, *68*, 108334.  
    https://doi.org/10.1016/j.rinp.2025.108334  
    Li2S gap ~3.2 eV in this GGA-level study; optical constants for the Li2X family. Use for qualitative trends (S → Se → Te: gap shrinks, reflectivity rises), not for quantitative T(λ).

Optical diagnostics to copy into any transparent-Li2S protocol:

- UV–vis T and R of the film on glass/ITO, 350–800 nm.
- Raman: T2g at ~373 cm−1 = crystalline Li2S; extra bands = polysulfide or amorphous chains (Klein).
- XPS S 2p: Li2S-type sulfide vs terminal/bridging polysulfide.
- Color: white = wide-gap Li2S; yellow = LiPS contamination.
- XRD: antifluorite, JCPDS 23-0369.

---

## E. Transparent current collectors, electrolytes, anodes

41. **Bates, J. B.; Dudney, N. J.; Gruzalski, G. R.; Zuhr, R. A.; Choudhury, A.; Luck, C. F.; Robertson, J. D.** Fabrication and characterization of amorphous lithium electrolyte thin films and rechargeable thin-film batteries. *J. Power Sources* **1993**, *43*, 103–110.  
    https://doi.org/10.1016/0378-7753(93)80106-Y  
    RF sputter of Li3PO4 in N2 → LiPON. σ ≈ 2×10−6 S cm−1 at 25 °C. Stable vs Li metal. The electrolyte of every inorganic transparent TFB since.

42. **Yu, X.; Bates, J. B.; Jellison, G. E.; Hart, F. X.** A stable thin-film lithium electrolyte: lithium phosphorus oxynitride. *J. Electrochem. Soc.* **1997**, *144*, 524–532.  
    https://doi.org/10.1149/1.1837443  
    Composition–conductivity map. Pair with Bates 1993.

43. **A review on lithium phosphorus oxynitride.** *J. Phys. Chem. C* **2021**.  
    https://doi.org/10.1021/acs.jpcc.0c10001  
    Stoichiometry, structure–conductivity, degradation vs Li, ALD vs sputter. Enough LiPON after Bates/Yu.

44. **Kozen, A. C.; Pearse, A. J.; Lin, C.-F.; Noked, M.; Rubloff, G. W.** Atomic layer deposition of the solid electrolyte LiPON. *Chem. Mater.* **2015**, *27*, 5324–5331.  
    https://doi.org/10.1021/acs.chemmater.5b01654  
    ALD LiPON for 3D / low-T / conformal coatings on transparent scaffolds.

45. **Lacivita, V. et al.** Resolving the amorphous structure of lithium phosphorus oxynitride (LiPON). *J. Am. Chem. Soc.* **2018**, *140*, 11029–11038.  
    https://doi.org/10.1021/jacs.8b05192  
    Amorphous structure that actually conducts.

46. **Xiao, Y.; Wang, Y.; Bo, S.-H.; Kim, J. C.; Miara, L.; Ceder, G.** Understanding interface stability in solid-state batteries. *Nat. Rev. Mater.* **2020**, *5*, 105–126.  
    https://doi.org/10.1038/s41578-019-0157-5  
    Thermodynamic interface map. Required before pairing Li2S with any SE other than LiPON.

47. **Manthiram, A.; Yu, X.; Wang, S.** Lithium battery chemistries enabled by solid-state electrolytes. *Nat. Rev. Mater.* **2017**, *2*, 16103.  
    https://doi.org/10.1038/natrevmats.2016.103  
    SE family survey (garnet, NASICON, sulfide, polymer, LiPON).

48. **Puthirath, A. B. et al.** Transparent flexible lithium ion conducting solid polymer electrolyte. *J. Mater. Chem. A* **2017**, *5*, 11152–11162.  
    https://doi.org/10.1039/C7TA02182H  
    Transparent SPE option if PVD LiPON is unavailable.

49. **Oukassi, S.; Giroud-Garampon, C.; Dubarry, C.; Ducros, C.; Salot, R.** All inorganic thin film electrochromic device using LiPON as the ion conductor. *Sol. Energy Mater. Sol. Cells* **2015**, *145*, 2–7.  
    https://doi.org/10.1016/j.solmat.2015.06.052  
    LiPON already proven in a transparent optical device (electrochromic), not only in batteries. Direct optical-stack precedent.

Anode constraint for a transparent Li2S full cell: Li metal is opaque and dendritic. Transparent options are (a) Li-free Si or IGZO thin film on TCO (Oukassi, Huang), (b) anode-less Cu/ITO with Li plated from Li2S, (c) lithiated transparent oxide. Li2S is the cathode that makes (b) and Li-free anodes chemically possible.

---

## F. Operando optical methods

Not device papers. They are how Li2S formation is watched through a window, and they set the spectroscopic background for a transparent cell that is also an analytical platform (Yang 2011 already noted in-situ Raman as a side benefit).

50. **In-operando imaging of polysulfide catholytes for Li–S batteries and implications for kinetics and mechanical stability.** *J. Power Sources* **2019**.  
    https://doi.org/10.1016/j.jpowsour.2019.05.038  
    Donut cell. Direct imaging of solid Li2S film growth from a liquid catholyte. Stress vs rate: high rate → less dense Li2S.

51. **Polysulfide speciation in Li–S battery electrolyte via in-operando optical imaging and ex-situ UV–vis spectra analysis.** *J. Electrochem. Soc.* **2022**.  
    https://doi.org/10.1149/1945-7111/ac8b3d  
    Color of the electrolyte as a LiPS-chain-length readout. Fully vs sparingly solvating electrolytes.

52. **Visualising the effect of areal current density on the performance and degradation of lithium sulfur batteries using operando optical microscopy.** *J. Electrochem. Soc.* **2024**.  
    https://doi.org/10.1149/1945-7111/ad9cc6  
    Sapphire-window cell. Current-density dependence of Li2S morphology.

53. **Direct tracking of the polysulfide shuttling and interfacial evolution in all-solid-state lithium–sulfur batteries.** *Energy Environ. Sci.* **2019**, *12*, 250–?  
    https://doi.org/10.1039/C9EE00578A  
    In-situ OM of polymer–ceramic ASSLS cells. Bright-white → dark-brown electrolyte = LiPS shuttle inside a “solid” electrolyte. Temperature dependence. Warning: a transparent solid-state Li2S cell with a polymer SE is not automatically shuttle-free.

54. **Beyene, T. T. et al.** Probing Li2S activation mechanism in lithium–sulfur batteries via multimodal operando techniques. *ACS Energy Lett.* **2026**.  
    https://doi.org/10.1021/acsenergylett.5c03358  
    Newest activation-mechanism paper. Read after Yang 2012 and Zhang 2017.

---

## G. Open data, structures, models

55. **Materials Project** Li2S entries mp-1153 (cubic) and mp-1125 (orthorhombic). CIF, XRD, band structure, elasticity.  
    https://next-gen.materialsproject.org/materials/mp-1153

56. **ComBat database** (Rashatwi et al.). Quantum-chemical and MD properties of Li–S electrolytes, including Li2S8 in DOL blends. PDB, binding energies, RDFs, diffusion.  
    https://github.com/rashatwi/combat

57. **Lithium-polysulfide structure set** (CIF, POSCAR, QE).  
    https://github.com/abineshperumal20-creator/Lithium-Polysulfides

58. **SPAN literature corpus** (sulfurized PAN cathode bibliography, 200+ entries). Adjacent chemistry, not Li2S, but the most complete open Li–S cathode citation dump.  
    https://github.com/weimufeng/SPAN

---

## H. Stack that a transparent Li2S cell would actually be

Inferred, not demonstrated:

```
glass or PET
  ITO or FTO or graphene / AgNW          # transparent current collector
  Li2S, 20–200 nm, sputtered or ALD      # Klein 2017 / Meng 2014
    optional: sub-eye-resolution carbon or VG grid if electronic wiring is required
  LiPON, 0.5–1.5 µm, RF sputter          # Bates 1993 / Oukassi 2019
  Si or IGZO or evaporated Li, patterned # Oukassi 2019 / Huang 2023 / Deng 2023
  ITO top contact, grid or thin
```

Optical budget: every opaque layer must be either (i) thinner than the absorption length in the visible or (ii) occupying <40% of the projected area as a grid. Li2S itself, if white/crystalline and thin, is not the optical bottleneck. Carbon hosts, Li metal, and thick ITO are.

Materials trade-off that will dominate the project: S vacancies and carbon raise electronic conductivity and destroy transmittance (Li 2017; Deng 2023 uses VGs). A transparent cell will likely need a lateral metal/TCO grid plus a thin, vacancy-poor Li2S film, accepting a first-charge activation protocol (Yang 2012) or a redox-mediating interlayer that is itself transparent.

---

## I. Full citation list (DOI order of appearance)

| ID | Citation | DOI |
|----|----------|-----|
| 1 | Yang et al., *PNAS* 2011 | 10.1073/pnas.1102873108 |
| 2 | Oukassi et al., *ACS AMI* 2019 | 10.1021/acsami.8b16364 |
| 3 | Yang et al., *JACS* 2012 | 10.1021/ja3052206 |
| 4 | Ting et al., *ACS Omega* 2022 | 10.1021/acsomega.2c05477 |
| 5 | Klein & Manthiram, *JACS* 2017 | 10.1021/jacs.7b03379 |
| 6 | Deng et al., *Nano-Micro Lett.* 2023 | 10.1007/s40820-023-01064-y |
| 7 | Bates et al., *J. Power Sources* 1993 | 10.1016/0378-7753(93)80106-Y |
| 8 | Eithiraj et al., *pss (b)* 2007 | 10.1002/pssb.200642506 |
| 9 | Wehner et al., *ACS Cent. Sci.* 2021 | 10.1021/acscentsci.0c01318 |
| 10 | Kim et al., *Nat. Chem. Eng.* 2024 | 10.1038/s44286-024-00079-5 |
| 11 | Pat et al., *Mater. Res. Express* 2019 | 10.1088/2053-1591/aae4aa |
| 12 | Huang et al., *Nat. Commun.* 2023 | 10.1038/s41467-023-41181-1 |
| 13 | Lee et al., *Sci. Rep.* 2022 | 10.1038/s41598-022-19933-8 |
| 14 | Kaiser et al., *Energy Storage Mater.* 2019 | 10.1016/j.ensm.2019.04.001 |
| 15 | Gao et al., *J. Mater. Chem. A* 2023 | 10.1039/D3TA06057H |
| 16 | Park et al., *Small* 2026 | 10.1002/smll.202513644 |
| 17 | Ji, Lee, Nazar, *Nat. Mater.* 2009 | 10.1038/nmat2460 |
| 18 | Manthiram et al., *Chem. Rev.* 2014 | 10.1021/cr500062v |
| 19 | Seh et al., *Chem. Soc. Rev.* 2016 | 10.1039/C5CS00410A |
| 20 | Pang et al., *Nat. Energy* 2016 | 10.1038/nenergy2016132 |
| 21 | Ould Ely et al., *ACS Appl. Energy Mater.* 2018 | 10.1021/acsaem.7b00153 |
| 21b | Zhou, Chen, Cui, *Nat. Energy* 2022 | 10.1038/s41560-022-01001-0 |
| 22 | Zhang et al., *Nano Lett.* 2017 | 10.1021/acs.nanolett.7b02381 |
| 23 | Zu, Klein, Manthiram, *J. Phys. Chem. Lett.* 2014 | 10.1021/jz5021108 |
| 24 | Ye et al., *ACS Energy Lett.* 2020 | 10.1021/acsenergylett.0c00936 |
| 25 | Tan et al., *Nat. Energy* 2017 | 10.1038/nenergy.2017.90 |
| 26 | Meng et al., *ACS Nano* 2014 | 10.1021/nn505480w |
| 27 | Wang, Jhang, Kou et al., *Nat. Commun.* 2023 | 10.1038/s41467-023-37564-z |
| 28 | Ekuma et al., *IJMPB* 2015 | 10.1142/S0217979215420060 |
| 29 | Bertheville et al., *J. Phys.: Condens. Matter* 1998 | 10.1088/0953-8984/10/9/018 |
| 30 | See et al., *J. Electrochem. Soc.* 2016 | 10.1149/2.0051610jes |
| 31 | Yu et al., *J. Electrochem. Soc.* 1997 | 10.1149/1.1837443 |
| 32 | Kozen et al., *Chem. Mater.* 2015 | 10.1021/acs.chemmater.5b01654 |
| 33 | Kozen/Put LiPON review, *J. Phys. Chem. C* 2021 | 10.1021/acs.jpcc.0c10001 |
| 34 | Lacivita et al., *JACS* 2018 | 10.1021/jacs.8b05192 |
| 35 | Xiao et al., *Nat. Rev. Mater.* 2020 | 10.1038/s41578-019-0157-5 |
| 36 | Manthiram, Yu, Wang, *Nat. Rev. Mater.* 2017 | 10.1038/natrevmats.2016.103 |
| 37 | Puthirath et al., *J. Mater. Chem. A* 2017 | 10.1039/C7TA02182H |
| 38 | Oukassi et al., *Sol. Energy Mater. Sol. Cells* 2015 | 10.1016/j.solmat.2015.06.052 |
| 39 | EEE 2019 in-operando imaging, *J. Power Sources* | 10.1016/j.jpowsour.2019.05.038 |
| 40 | *JES* 2022 optical + UV–vis LiPS | 10.1149/1945-7111/ac8b3d |
| 41 | *JES* 2024 operando OM vs current density | 10.1149/1945-7111/ad9cc6 |
| 42 | *Energy Environ. Sci.* 2019 ASSLS optical shuttle | 10.1039/C9EE00578A |
| 43 | Nagao, Hayashi, Tatsumisago, *J. Mater. Chem.* 2012 | 10.1039/c2jm16802b |
| 44 | Li vacancy-gap paper, *Funct. Mater. Lett.* 2017 | 10.1142/S1793604717500679 |
| 45 | *Results Phys.* 2025 Li2X optical DFT | 10.1016/j.rinp.2025.108334 |
| 46 | Su et al., *Adv. Sustainable Syst.* 2025 | 10.1002/adsu.202400555 |
| 47 | Liang et al., *Adv. Energy Mater.* 2024 | 10.1002/aenm.202401959 |
| 48 | *Nano Res. Energy* 2025 Li2S host/electrolyte | 10.26599/NRE.2025.9120194 |
| 49 | Beyene et al., *ACS Energy Lett.* 2026 | 10.1021/acsenergylett.5c03358 |
| 50 | Klein dissertation, UT Austin 2016 | 10.15781/t2x63bb46 |

Open-access starting set (no paywall): Yang 2011 (PMC), Ting 2022 (PMC), Deng 2023 (Springer OA), Park 2026 (Chalmers PDF), Yang 2012 (Cui group PDF), Klein dissertation (UT Austin), Bertheville 1998 (Geneva PDF), Wehner 2021 (PMC), Jia 2023 (PMC).

---

## J. What not to treat as a transparent Li2S battery paper

- Operando optical Li–S cells with sapphire windows. They are analytical cells, not transparent devices.
- Bulk all-solid-state Li2S pellets under 50–300 MPa. Opaque by construction.
- Carbon-hosted Li2S nanocomposites (graphene, CNT, KB). Optically black.
- Commercial Li2S powder papers that do not report UV–vis. Yellow powder is not a transparent film.

A paper counts as on-topic for *transparent Li2S* only if it reports (i) visible transmittance of a Li2S-containing stack, or (ii) a vapor-phase Li2S film whose thickness and substrate are compatible with T > 40%, or (iii) an architecture that has already been shown transparent in Li-ion chemistry and is chemically substitutable with Li2S. Category (i) is currently empty.
