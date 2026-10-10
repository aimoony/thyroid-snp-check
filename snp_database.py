"""
Curated SNP database for snp_report.py.

All alleles are on the GRCh37 PLUS (forward) strand, which is how 23andMe,
AncestryDNA (mostly) and most consumer raw-data files report genotypes.

Fields
------
rsid      : identifier as it appears in the raw file (rs… or 23andMe internal i…)
gene      : gene symbol
variant   : common name of the variant
cat       : report section key (see CATEGORIES)
alleles   : the two plus-strand alleles expected at this position
effect    : the allele the description refers to (counted 0/1/2 in your genotype)
kind      : "risk"    – effect allele is the less favourable / variant allele
            "benefit" – effect allele is the favourable allele
            "info"    – direction unclear or purely descriptive
evidence  : "strong" | "moderate" | "limited"
desc      : what carrying the effect allele is associated with
notes     : optional {copies: text} overrides for 0 / 1 / 2 copies
"""

CATEGORIES = [
    # key, title, group
    ("thy_conv", "T4 → T3 conversion (deiodinases)", "Thyroid"),
    ("thy_tsh", "TSH set-point & thyroid function", "Thyroid"),
    ("thy_auto", "Thyroid autoimmunity (Hashimoto's / Graves')", "Thyroid"),
    ("thy_cancer", "Thyroid nodule / cancer susceptibility", "Thyroid"),
    ("thy_cof", "Thyroid cofactors (selenium, iron, vitamin A, D, B12)", "Thyroid"),
    ("methyl", "Methylation & folate cycle", "Common"),
    ("cardio", "Cardiovascular, lipids & clotting", "Common"),
    ("metab", "Metabolic, weight & liver", "Common"),
    ("pgx", "Drug response (pharmacogenomics)", "Common"),
    ("detox", "Detox & oxidative stress", "Common"),
    ("nutri", "Nutrition, neuro & traits", "Common"),
    ("carrier", "Other health / carrier variants", "Common"),
]


def S(rsid, gene, variant, cat, alleles, effect, kind, evidence, desc, notes=None):
    return dict(rsid=rsid, gene=gene, variant=variant, cat=cat, alleles=alleles,
                effect=effect, kind=kind, evidence=evidence, desc=desc,
                notes=notes or {})


SNPS = [
    # ------------------------------------------------------------------ THYROID
    # --- conversion
    S("rs225014", "DIO2", "Thr92Ala", "thy_conv", "CT", "C", "risk", "moderate",
      "Ala92 allele. Linked to modestly reduced type-2 deiodinase activity (local T4→T3 conversion "
      "in brain, muscle, pituitary). In levothyroxine-treated patients, Ala/Ala has been associated with "
      "lower well-being and a better response to T4+T3 combination therapy in some studies; results "
      "are mixed and current guidelines do not yet use it to choose therapy.",
      {2: "Ala/Ala — the genotype most often discussed for possible benefit from T4+T3 combination therapy.",
       1: "Thr/Ala — one copy; intermediate, generally small effect.",
       0: "Thr/Thr — common, typical D2 activity."}),
    S("rs12885300", "DIO2", "ORFa-Gly3Asp", "thy_conv", "CT", "T", "info", "limited",
      "Asp3 allele of an upstream open reading frame in DIO2. Associated with slightly lower T4, FT4 and "
      "rT3 (suggesting somewhat higher D2 activity; Peeters 2005); effects are small and not consistently replicated."),
    S("rs2235544", "DIO1", "rs2235544 A>C", "thy_conv", "AC", "C", "benefit", "moderate",
      "C allele associated with higher type-1 deiodinase activity: higher FT3/FT4 ratio, lower FT4 and "
      "lower reverse T3 (Panicker 2008). The A allele is linked to relatively less T4→T3 conversion by D1.",
      {0: "A/A — relatively lower D1 conversion (higher FT4, lower FT3/FT4 ratio on average).",
       1: "A/C — intermediate.",
       2: "C/C — highest average D1 activity (higher FT3/FT4 ratio)."}),
    S("rs11206244", "DIO1", "D1-C785T", "thy_conv", "CT", "T", "risk", "moderate",
      "T allele associated with higher reverse T3 and a lower T3/rT3 and FT3/FT4 ratio, i.e. slightly "
      "less efficient D1 conversion (Peeters 2003; Panicker 2008)."),
    S("rs12095080", "DIO1", "D1a-A/G", "thy_conv", "AG", "G", "benefit", "moderate",
      "G allele associated with a higher T3/rT3 ratio (more active D1) (Peeters 2003).",
      {0: "A/A — common genotype; typical D1 activity.",
       1: "A/G — one copy of the higher-T3/rT3 allele.",
       2: "G/G — two copies of the higher-T3/rT3 allele."}),

    # --- TSH / thyroid function
    S("rs4704397", "PDE8B", "rs4704397", "thy_tsh", "AG", "A", "info", "strong",
      "A allele associated with higher serum TSH (≈ +0.1 mU/L per copy; Arnaud-Lopez 2008). Affects the "
      "TSH reference point, not thyroid disease per se — a mildly higher TSH can be genetic."),
    S("rs6885099", "PDE8B", "rs6885099", "thy_tsh", "AG", "G", "info", "strong",
      "G allele associated with higher serum TSH (A lowers it; Porcu 2013) and modestly with hypothyroidism."),
    S("rs7850258", "FOXE1", "rs7850258", "thy_tsh", "AG", "G", "risk", "moderate",
      "G allele associated with higher risk of primary hypothyroidism (Denny 2011 PheWAS; Eriksson 2012 "
      "23andMe GWAS). The minor A allele is protective (OR ≈ 0.75)."),
    S("rs1991517", "TSHR", "Asp727Glu", "thy_tsh", "CG", "G", "info", "limited",
      "Glu727 allele; associated with slightly lower TSH in some cohorts. Effect is small."),
    S("rs4915077", "VAV3", "rs4915077", "thy_tsh", "CT", "C", "risk", "moderate",
      "Minor C allele (~7% in Europeans) associated with hypothyroidism (OR ≈ 1.3; Eriksson 2012 23andMe GWAS)."),

    # --- autoimmunity
    S("rs2476601", "PTPN22", "R620W (1858C>T)", "thy_auto", "AG", "A", "risk", "strong",
      "620W allele; one of the best-replicated autoimmunity variants (Hashimoto's, Graves', type 1 "
      "diabetes, rheumatoid arthritis). OR ≈ 1.5–2 per copy for autoimmune thyroid disease."),
    S("rs6679677", "PTPN22 region", "rs6679677", "thy_auto", "AC", "A", "risk", "strong",
      "Proxy for PTPN22 R620W; A allele associated with hypothyroidism in the 23andMe GWAS."),
    S("rs3087243", "CTLA4", "CT60", "thy_auto", "AG", "G", "risk", "moderate",
      "G allele associated with Graves' disease, Hashimoto's and thyroid antibodies (Ueda 2003). "
      "G is the more common allele in most populations, so G/G is frequent.",
      ),
    S("rs231775", "CTLA4", "+49A/G (Thr17Ala)", "thy_auto", "AG", "G", "risk", "moderate",
      "Ala17 (G) allele associated with Graves' disease and autoimmune hypothyroidism (OR ≈ 1.3–1.5)."),
    S("rs3184504", "SH2B3", "R262W", "thy_auto", "CT", "T", "risk", "strong",
      "T allele associated with hypothyroidism, TPO antibodies, celiac disease and type 1 diabetes."),
    S("rs7574865", "STAT4", "rs7574865", "thy_auto", "GT", "T", "risk", "moderate",
      "T allele associated with several autoimmune diseases (lupus, RA) and with autoimmune thyroid disease."),
    S("rs7528684", "FCRL3", "-169T>C", "thy_auto", "AG", "G", "risk", "moderate",
      "Risk allele associated with Graves' disease and RA in Japanese and some European cohorts."),
    S("rs2187668", "HLA-DQA1", "DQ2.5 tag", "thy_auto", "CT", "T", "risk", "strong",
      "T allele tags HLA-DQ2.5, the main genetic risk factor for celiac disease. Celiac and Hashimoto's "
      "frequently co-occur."),
    S("rs7454108", "HLA-DQB1", "DQ8 tag", "thy_auto", "CT", "C", "risk", "strong",
      "C allele tags HLA-DQ8, the second celiac-risk haplotype; also linked to type 1 diabetes and "
      "autoimmune thyroid disease. Carrying DQ8 is common (~10–20%) and most carriers never develop celiac."),
    S("rs179247", "TSHR", "intron 1", "thy_auto", "AG", "A", "risk", "strong",
      "A allele associated with Graves' disease (TSHR intron-1 locus; Dechairo 2005)."),
    S("rs12101255", "TSHR", "intron 1", "thy_auto", "CT", "T", "risk", "strong",
      "T allele associated with Graves' disease (Brand 2009)."),
    S("rs2268458", "TSHR", "intron 1", "thy_auto", "CT", "C", "risk", "limited",
      "C allele associated with Graves' disease at the TSHR intron-1 locus (single small study; Yin 2008)."),

    # --- cancer / nodules
    S("rs965513", "FOXE1/PTCSC2", "rs965513", "thy_cancer", "AG", "A", "risk", "strong",
      "A allele associated with papillary thyroid cancer (OR ≈ 1.7 per copy; Gudmundsson 2009) and "
      "with lower TSH / T4. Absolute risk remains low."),
    S("rs1867277", "FOXE1", "promoter -283G>A", "thy_cancer", "AG", "A", "risk", "moderate",
      "A allele associated with papillary thyroid cancer (Landa 2009)."),
    S("rs944289", "PTCSC3 / NKX2-1", "rs944289", "thy_cancer", "CT", "T", "risk", "strong",
      "T allele associated with papillary thyroid cancer (OR ≈ 1.4; Gudmundsson 2009)."),
    S("rs966423", "DIRC3", "rs966423", "thy_cancer", "CT", "C", "risk", "moderate",
      "C allele associated with differentiated thyroid cancer and with TSH levels (Gudmundsson 2012)."),

    # --- cofactors
    S("rs1050450", "GPX1", "Pro198Leu", "thy_cof", "AG", "A", "risk", "moderate",
      "Leu198 allele; lower glutathione peroxidase-1 activity (a selenium-dependent enzyme). The thyroid "
      "relies on selenoproteins to handle the H2O2 it uses to make hormone."),
    S("rs3877899", "SELENOP", "Ala234Thr", "thy_cof", "CT", "T", "info", "limited",
      "Selenoprotein P variant; alters plasma selenium / selenoprotein P levels and response to "
      "selenium supplementation."),
    S("rs7579", "SELENOP", "3'UTR 25191G>A", "thy_cof", "CT", "T", "info", "limited",
      "Selenoprotein P 3'UTR variant; associated with differences in selenoprotein P levels and "
      "response to selenium supplementation."),
    S("rs855791", "TMPRSS6", "Ala736Val", "thy_cof", "AG", "A", "risk", "strong",
      "Val736 allele associated with lower serum iron, transferrin saturation and haemoglobin. Iron is "
      "required by thyroid peroxidase; low ferritin is common in hypothyroidism."),
    S("rs1800562", "HFE", "C282Y", "thy_cof", "AG", "A", "risk", "strong",
      "Main hereditary haemochromatosis (iron overload) variant. Two copies = high risk of iron overload "
      "(check ferritin / transferrin saturation); one copy = carrier."),
    S("rs1799945", "HFE", "H63D", "thy_cof", "CG", "G", "risk", "moderate",
      "Milder haemochromatosis variant; matters mainly together with C282Y."),
    S("rs12934922", "BCO1 (BCMO1)", "R267S", "thy_cof", "AT", "T", "risk", "moderate",
      "Reduced conversion of beta-carotene to vitamin A (Leung 2009). Vitamin A influences TSH and "
      "thyroid-hormone receptor signalling."),
    S("rs7501331", "BCO1 (BCMO1)", "A379V", "thy_cof", "CT", "T", "risk", "moderate",
      "Reduced beta-carotene → vitamin A conversion; combined with R267S can lower conversion ~50–70%."),
    S("rs2282679", "GC (vitamin D binding protein)", "rs2282679", "thy_cof", "GT", "G", "risk", "strong",
      "Allele associated with lower circulating 25-OH vitamin D (Wang 2010). Low vitamin D is "
      "associated with thyroid autoimmunity."),
    S("rs1544410", "VDR", "BsmI", "thy_cof", "CT", "T", "info", "limited",
      "T = 'B' allele (in strong LD with TaqI 't'). Associations with bone density and autoimmunity are "
      "weak/inconsistent."),
    S("rs731236", "VDR", "TaqI", "thy_cof", "AG", "G", "info", "limited",
      "G = 't' allele. Associations with vitamin-D response and autoimmune thyroid disease are weak/inconsistent."),
    S("rs601338", "FUT2", "W143X (non-secretor)", "thy_cof", "AG", "A", "info", "strong",
      "A/A = 'non-secretor': higher serum B12, altered gut microbiome, resistance to norovirus."),
    S("rs1801198", "TCN2", "C776G (Pro259Arg)", "thy_cof", "CG", "G", "risk", "moderate",
      "776G (Arg259) allele associated with lower holotranscobalamin (active B12 delivered to cells). "
      "B12 deficiency is common in autoimmune hypothyroidism. (One large proteomics study reported the "
      "opposite direction, likely an assay artefact.)"),

    # ------------------------------------------------------------------ METHYLATION
    S("rs1801133", "MTHFR", "C677T (Ala222Val)", "methyl", "AG", "A", "risk", "strong",
      "677T allele. T/T reduces MTHFR activity to ~30% and can raise homocysteine, especially with low "
      "folate/riboflavin; C/T ~65% activity."),
    S("rs1801131", "MTHFR", "A1298C (Glu429Ala)", "methyl", "GT", "G", "risk", "moderate",
      "1298C allele. Milder than C677T: C/C reduces activity modestly and usually does not raise "
      "homocysteine on its own."),
    S("rs1805087", "MTR", "A2756G (D919G)", "methyl", "AG", "G", "info", "limited",
      "2756G allele of methionine synthase (B12-dependent). Effects on homocysteine are small."),
    S("rs1801394", "MTRR", "A66G (I22M)", "methyl", "AG", "G", "risk", "limited",
      "66G allele; slightly reduced methionine synthase reductase efficiency (B12 recycling)."),
    S("rs2236225", "MTHFD1", "G1958A (R653Q)", "methyl", "AG", "A", "risk", "limited",
      "653Q allele; associated with greater choline/folate demand and neural-tube-defect risk in some studies."),
    S("rs234706", "CBS", "C699T (Y233Y)", "methyl", "AG", "A", "info", "limited",
      "699T allele. Popularly called 'CBS upregulation', but the evidence is weak; associated with "
      "slightly lower homocysteine."),
    S("rs1979277", "SHMT1", "C1420T (L474F)", "methyl", "AG", "A", "info", "limited",
      "1420T allele; small effects on folate partitioning."),
    S("rs3733890", "BHMT", "G742A (R239Q)", "methyl", "AG", "A", "info", "limited",
      "239Q allele of betaine-homocysteine methyltransferase (alternate, choline/betaine-based "
      "homocysteine pathway)."),
    S("rs7946", "PEMT", "V175M", "methyl", "CT", "T", "risk", "limited",
      "Variant associated with reduced PEMT activity and higher dietary choline requirement."),
    S("rs4680", "COMT", "Val158Met", "methyl", "AG", "A", "info", "strong",
      "Met158 allele ≈ 3–4× lower COMT activity (slower breakdown of dopamine, adrenaline, catechol-estrogens).",
      {0: "Val/Val — 'fast' COMT: quicker dopamine/catecholamine clearance; often more stress-resilient, "
          "sometimes lower baseline prefrontal dopamine.",
       1: "Val/Met — intermediate COMT activity.",
       2: "Met/Met — 'slow' COMT: higher prefrontal dopamine, slower estrogen/catecholamine clearance."}),
    S("rs6323", "MAOA", "R297R (T941G)", "methyl", "GT", "T", "info", "limited",
      "Synonymous MAOA variant on the X chromosome (men have one copy). T is usually reported as the "
      "lower-MAOA-activity allele and G as higher, but the literature is inconsistent; clinical meaning is uncertain."),

    # ------------------------------------------------------------------ CARDIO
    S("rs6025", "F5", "Factor V Leiden (R506Q)", "cardio", "CT", "T", "risk", "strong",
      "Factor V Leiden. One copy ≈ 3–8× venous clot risk; two copies much higher."),
    S("i3002432", "F2", "Prothrombin G20210A", "cardio", "AG", "A", "risk", "strong",
      "Prothrombin G20210A (23andMe internal ID). One copy ≈ 2–3× venous clot risk."),
    S("rs10757274", "CDKN2B-AS1 (9p21)", "9p21 CAD locus", "cardio", "AG", "G", "risk", "strong",
      "G allele associated with coronary artery disease (~1.2–1.3× per copy)."),
    S("rs1333049", "CDKN2B-AS1 (9p21)", "9p21 CAD locus", "cardio", "CG", "C", "risk", "strong",
      "C allele associated with coronary artery disease (~1.2–1.3× per copy)."),
    S("rs10455872", "LPA", "Lp(a) variant", "cardio", "AG", "G", "risk", "strong",
      "G allele strongly raises lipoprotein(a); ~1.5× coronary risk per copy."),
    S("rs3798220", "LPA", "I4399M", "cardio", "CT", "C", "risk", "strong",
      "C allele raises lipoprotein(a) and coronary risk."),
    S("rs5742904", "APOB", "R3527Q (familial defective ApoB)", "cardio", "CT", "T", "risk", "strong",
      "Causes familial defective ApoB-100 (high LDL)."),
    S("rs11591147", "PCSK9", "R46L", "cardio", "GT", "T", "benefit", "strong",
      "Loss-of-function allele: lifelong lower LDL cholesterol and ~45% lower coronary risk."),
    S("rs4343", "ACE", "I/D proxy (2350G>A)", "cardio", "AG", "G", "info", "limited",
      "G tags the ACE 'D' allele (higher ACE levels). Associations with blood pressure/athletic power are modest."),

    # ------------------------------------------------------------------ METABOLIC
    S("rs7903146", "TCF7L2", "rs7903146", "metab", "CT", "T", "risk", "strong",
      "Strongest common type 2 diabetes variant (~1.4× per copy)."),
    S("rs10811661", "CDKN2A/B", "rs10811661", "metab", "CT", "T", "risk", "strong",
      "T allele associated with type 2 diabetes (~1.2×)."),
    S("rs13266634", "SLC30A8", "R325W", "metab", "CT", "C", "risk", "strong",
      "C allele (Arg325) associated with type 2 diabetes (~1.1–1.2×)."),
    S("rs10830963", "MTNR1B", "rs10830963", "metab", "CG", "G", "risk", "strong",
      "G allele associated with higher fasting glucose and T2D risk; late-night eating may matter more for carriers."),
    S("rs1801282", "PPARG", "Pro12Ala", "metab", "CG", "G", "benefit", "strong",
      "Ala12 allele associated with modestly lower type 2 diabetes risk and better insulin sensitivity."),
    S("rs9939609", "FTO", "rs9939609", "metab", "AT", "A", "risk", "strong",
      "A allele: ~+0.4 kg/m² BMI per copy on average and higher appetite; effect is blunted by physical activity."),
    S("rs17782313", "MC4R", "rs17782313", "metab", "CT", "C", "risk", "strong",
      "C allele associated with higher BMI and appetite."),
    S("rs738409", "PNPLA3", "I148M", "metab", "CG", "G", "risk", "strong",
      "148M allele; the strongest genetic risk factor for fatty liver disease (NAFLD/MASLD) and higher ALT. "
      "Risk is amplified by excess weight, alcohol and high sugar/fructose intake.",
      {2: "M/M — roughly 2–3× higher risk of fatty liver and fibrosis; worth knowing for liver-enzyme monitoring.",
       1: "I/M — one copy; moderately increased fatty-liver susceptibility.",
       0: "I/I — typical risk."}),

    # ------------------------------------------------------------------ PGx
    S("rs4244285", "CYP2C19", "*2", "pgx", "AG", "A", "risk", "strong",
      "Loss-of-function. Reduces activation of clopidogrel (Plavix); alters dosing of many SSRIs and PPIs."),
    S("rs4986893", "CYP2C19", "*3", "pgx", "AG", "A", "risk", "strong",
      "Loss-of-function (mainly East Asian)."),
    S("rs12248560", "CYP2C19", "*17", "pgx", "CT", "T", "info", "strong",
      "Increased-function allele: faster metabolism of PPIs, some SSRIs (e.g. escitalopram, sertraline)."),
    S("rs1799853", "CYP2C9", "*2", "pgx", "CT", "T", "risk", "strong",
      "Reduced function: lower warfarin and NSAID clearance."),
    S("rs1057910", "CYP2C9", "*3", "pgx", "AC", "C", "risk", "strong",
      "Strongly reduced function: lower warfarin dose needed; slower NSAID / phenytoin clearance."),
    S("rs9923231", "VKORC1", "-1639G>A", "pgx", "CT", "T", "risk", "strong",
      "Increased warfarin sensitivity (lower dose needed)."),
    S("rs4149056", "SLCO1B1", "*5 (V174A)", "pgx", "CT", "C", "risk", "strong",
      "Higher simvastatin blood levels and statin-myopathy risk (C/C ≈ 17×)."),
    S("rs3892097", "CYP2D6", "*4", "pgx", "CT", "T", "risk", "strong",
      "Most common CYP2D6 non-functional allele (codeine, tramadol, tamoxifen, many antidepressants). "
      "Note: CYP2D6 copy-number changes cannot be seen on this chip."),
    S("rs3918290", "DPYD", "*2A", "pgx", "CT", "T", "risk", "strong",
      "Severe 5-FU / capecitabine toxicity risk."),
    S("rs1142345", "TPMT", "*3C (A719G)", "pgx", "CT", "C", "risk", "strong",
      "Reduced TPMT activity: azathioprine / 6-MP toxicity risk."),
    S("rs1800460", "TPMT", "*3B (G460A)", "pgx", "CT", "T", "risk", "strong",
      "Reduced TPMT activity: azathioprine / 6-MP toxicity risk."),
    S("rs776746", "CYP3A5", "*3", "pgx", "CT", "C", "info", "strong",
      "*3 = non-expresser. *3/*3 is the most common genotype in Europeans; relevant to tacrolimus dosing."),
    S("rs762551", "CYP1A2", "*1F (-163C>A)", "pgx", "AC", "C", "risk", "moderate",
      "C allele = slower inducible CYP1A2 (caffeine) metabolism. In slow metabolizers, heavy coffee intake "
      "has been linked to higher blood pressure / heart-attack risk.",
      {0: "A/A — fast caffeine metabolizer.",
       1: "A/C — slow(er) caffeine metabolizer.",
       2: "C/C — slow caffeine metabolizer."}),
    S("rs1800497", "ANKK1/DRD2", "Taq1A", "pgx", "AG", "A", "info", "moderate",
      "A1 allele associated with lower striatal D2 receptor density; reward/addiction-related research associations."),
    S("rs1799971", "OPRM1", "A118G (N40D)", "pgx", "AG", "G", "info", "moderate",
      "Asp40 allele; associated with higher opioid dose requirements and altered response to naltrexone."),

    # ------------------------------------------------------------------ DETOX
    S("rs1695", "GSTP1", "Ile105Val", "detox", "AG", "G", "info", "limited",
      "Val105 allele alters glutathione S-transferase P1 substrate activity."),
    S("rs4880", "SOD2", "Val16Ala", "detox", "AG", "G", "info", "limited",
      "Ala16 allele: more efficient import of MnSOD into mitochondria. Associations with disease are mixed."),
    S("rs1800566", "NQO1", "P187S", "detox", "AG", "A", "risk", "moderate",
      "Ser187 allele; T/T (here A/A) has almost no NQO1 activity (quinone / benzene detox, CoQ10 reduction)."),
    S("rs1800629", "TNF", "-308G>A", "detox", "AG", "A", "risk", "limited",
      "A allele associated with higher TNF-α production (inflammation)."),

    # ------------------------------------------------------------------ NUTRITION / TRAITS
    S("rs4988235", "MCM6/LCT", "-13910C>T", "nutri", "AG", "A", "info", "strong",
      "A allele (= T on the gene strand) gives lactase persistence (digesting lactose as an adult).",
      {0: "G/G — likely lactase non-persistent (lactose intolerance tendency in adulthood; "
          "note other persistence variants exist in African / Middle-Eastern populations).",
       1: "A/G — lactase persistent.",
       2: "A/A — lactase persistent."}),
    S("rs671", "ALDH2", "*2 (E504K)", "nutri", "AG", "A", "risk", "strong",
      "Alcohol flush; carriers who drink have much higher oesophageal-cancer risk."),
    S("rs1815739", "ACTN3", "R577X", "nutri", "CT", "T", "info", "moderate",
      "X allele (stop codon). X/X lacks alpha-actinin-3 in fast-twitch muscle (less sprint/power, possibly "
      "more endurance-oriented)."),
    S("rs174537", "FADS1", "rs174537", "nutri", "GT", "T", "risk", "strong",
      "T allele associated with lower conversion of plant omega-3/6 into EPA/DHA and arachidonic acid."),
    S("rs6265", "BDNF", "Val66Met", "nutri", "CT", "T", "info", "moderate",
      "Met66 allele: reduced activity-dependent BDNF secretion; memory and mood research associations."),
    S("rs12913832", "HERC2/OCA2", "eye colour", "nutri", "AG", "G", "info", "strong",
      "G/G strongly predicts blue eyes; A/A brown; A/G usually brown/hazel."),
    S("rs17822931", "ABCC11", "G538A", "nutri", "CT", "T", "info", "strong",
      "T/T = dry earwax and low underarm odour (common in East Asians)."),
    S("rs16969968", "CHRNA5", "D398N", "nutri", "AG", "A", "risk", "strong",
      "A allele associated with heavier smoking and nicotine dependence."),

    # ------------------------------------------------------------------ CARRIER / OTHER
    S("rs429358", "APOE", "C130R (ε4 defining SNP)", "carrier", "CT", "C", "risk", "strong",
      "See the APOE result in the summary — APOE type is determined from rs429358 + rs7412 together.",
      {1: "One ε4-defining allele (see APOE type in the summary).",
       2: "Two ε4-defining alleles (see APOE type in the summary)."}),
    S("rs7412", "APOE", "R176C (ε2 defining SNP)", "carrier", "CT", "T", "info", "strong",
      "See the APOE result in the summary."),
    S("rs28929474", "SERPINA1", "PiZ (E342K)", "carrier", "CT", "T", "risk", "strong",
      "Alpha-1 antitrypsin Z allele. Z/Z causes AAT deficiency (lung and liver disease)."),
    S("rs17580", "SERPINA1", "PiS (E264V)", "carrier", "AT", "A", "risk", "strong",
      "Alpha-1 antitrypsin S allele (mild); significant mainly as S/Z."),
]

# APOE haplotype lookup: (rs429358 genotype sorted, rs7412 genotype sorted) -> type
APOE_TABLE = {
    ("TT", "TT"): "ε2/ε2",
    ("TT", "CT"): "ε2/ε3",
    ("CT", "CT"): "ε2/ε4 (or ε1/ε3 — rare)",
    ("TT", "CC"): "ε3/ε3",
    ("CT", "CC"): "ε3/ε4",
    ("CC", "CC"): "ε4/ε4",
    ("CC", "CT"): "ε4/ε4 with ε1 (rare)",
}

APOE_NOTES = {
    "ε2/ε2": "Lower Alzheimer's risk; associated with type III hyperlipoproteinaemia in some carriers.",
    "ε2/ε3": "Slightly lower-than-average Alzheimer's risk; LDL often a bit lower.",
    "ε3/ε3": "The most common APOE type (≈ 60%); average risk.",
    "ε3/ε4": "One ε4 copy (≈ 25% of people). Associated with higher LDL cholesterol and roughly 2–3× "
             "average late-onset Alzheimer's risk — a risk factor, not a diagnosis; most ε3/ε4 carriers do "
             "not develop Alzheimer's. Cardiovascular fitness, blood-pressure, LDL, sleep and glucose "
             "control are the modifiable factors with the best evidence.",
    "ε4/ε4": "Two ε4 copies (≈ 2%). Substantially higher late-onset Alzheimer's risk; consider genetic counselling.",
}

# Genes whose full regions are exported to thyroid_gene_region_snps.csv
# (every SNP on the chip that falls inside these genes ± flank), so you get
# *all* measured thyroid-gene SNPs, not only the curated, annotated ones.
THYROID_GENES = [
    "DIO1", "DIO2", "DIO3", "TSHR", "TSHB", "TRH", "TRHR", "TPO", "TG", "IYD",
    "DUOX1", "DUOX2", "DUOXA2", "SLC5A5", "SLC26A4", "SLC26A7", "SLC16A2",
    "SLC16A10", "SLCO1C1", "THRA", "THRB", "PAX8", "NKX2-1", "FOXE1", "HHEX",
    "SECISBP2", "SEPP1", "GPX1", "GPX3", "TXNRD1", "PDE8B", "PDE10A",
    "CAPZB", "CTLA4", "PTPN22", "SH2B3", "VAV3", "FCRL3", "CD40", "IL2RA",
    "TTR", "SERPINA7", "ALB", "CRYM",
]
