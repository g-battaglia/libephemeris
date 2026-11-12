//! Body identifiers, physical constants, flags and public result slots.
//!
//! The project constants catalogue originates from LibEphemeris v3.2.1.
//! Integer selectors and flag combinations are exact; floating-point values
//! are stored as explicit binary64 representations. Table values are data,
//! not alternative implementations of the calculations that consume them.
//!
//! Identifier families cover planets, fixed stars, asteroids, hypothetical
//! bodies and planetary moons. Other tables define NAIF/Horizons mappings,
//! calendar and sidereal selectors, eclipse flags, and result indices.
//!
//! Units are documented per constant: epochs and spans use days, angular
//! quantities use their stated degrees/radians convention, distances use AU
//! or the indicated SI unit, and flags/identifiers are dimensionless.
//! Lookup helpers return `None` or `false` for an unmapped identifier;
//! physical body-domain validation belongs to the consuming calculation.

// Parentheses make identifier offsets and composite flag definitions
// explicit. The catalogue preserves public identifier values.
#![allow(unused_parens)]
#![allow(clippy::identity_op)]

/// Number of scalar/string constants in the base catalogue.
/// The four lookup functions are counted separately.
#[allow(dead_code)]
pub const RELEASE_CONST_COUNT: usize = 594;

/// Release constant `ECL_NUT`.
pub const ECL_NUT: i32 = -1;
/// Release constant `SUN`.
pub const SUN: i32 = 0;
/// Release constant `MOON`.
pub const MOON: i32 = 1;
/// Release constant `MERCURY`.
pub const MERCURY: i32 = 2;
/// Release constant `VENUS`.
pub const VENUS: i32 = 3;
/// Release constant `MARS`.
pub const MARS: i32 = 4;
/// Release constant `JUPITER`.
pub const JUPITER: i32 = 5;
/// Release constant `SATURN`.
pub const SATURN: i32 = 6;
/// Release constant `URANUS`.
pub const URANUS: i32 = 7;
/// Release constant `NEPTUNE`.
pub const NEPTUNE: i32 = 8;
/// Release constant `PLUTO`.
pub const PLUTO: i32 = 9;
/// Release constant `MEAN_NODE`.
pub const MEAN_NODE: i32 = 10;
/// Release constant `TRUE_NODE`.
pub const TRUE_NODE: i32 = 11;
/// Release constant `MEAN_APOG`.
pub const MEAN_APOG: i32 = 12;
/// Release constant `OSCU_APOG`.
pub const OSCU_APOG: i32 = 13;
/// Release constant `EARTH`.
pub const EARTH: i32 = 14;
/// Release constant `CHIRON`.
pub const CHIRON: i32 = 15;
/// Release constant `PHOLUS`.
pub const PHOLUS: i32 = 16;
/// Release constant `CERES`.
pub const CERES: i32 = 17;
/// Release constant `PALLAS`.
pub const PALLAS: i32 = 18;
/// Release constant `JUNO`.
pub const JUNO: i32 = 19;
/// Release constant `VESTA`.
pub const VESTA: i32 = 20;
/// Release constant `INTP_APOG`.
pub const INTP_APOG: i32 = 21;
/// Release constant `INTP_PERG`.
pub const INTP_PERG: i32 = 22;
/// Release constant `NPLANETS`.
pub const NPLANETS: i32 = 23;
/// Release constant `AST_OFFSET`.
pub const AST_OFFSET: i32 = 10000;
/// Release constant `VARUNA`.
pub const VARUNA: i32 = (AST_OFFSET + 20000);
/// Release constant `FICT_OFFSET`.
pub const FICT_OFFSET: i32 = 40;
/// Release constant `NFICT_ELEM`.
pub const NFICT_ELEM: i32 = 15;
/// Release constant `COMET_OFFSET`.
pub const COMET_OFFSET: i32 = 1000;
/// Release constant `CUPIDO`.
pub const CUPIDO: i32 = (FICT_OFFSET + 0);
/// Release constant `HADES`.
pub const HADES: i32 = (FICT_OFFSET + 1);
/// Release constant `ZEUS`.
pub const ZEUS: i32 = (FICT_OFFSET + 2);
/// Release constant `KRONOS`.
pub const KRONOS: i32 = (FICT_OFFSET + 3);
/// Release constant `APOLLON`.
pub const APOLLON: i32 = (FICT_OFFSET + 4);
/// Release constant `ADMETOS`.
pub const ADMETOS: i32 = (FICT_OFFSET + 5);
/// Release constant `VULKANUS`.
pub const VULKANUS: i32 = (FICT_OFFSET + 6);
/// Release constant `POSEIDON`.
pub const POSEIDON: i32 = (FICT_OFFSET + 7);
/// Release constant `ISIS`.
pub const ISIS: i32 = (FICT_OFFSET + 8);
/// Release constant `TRANSPLUTO`.
pub const TRANSPLUTO: i32 = ISIS;
/// Release constant `VULCAN`.
pub const VULCAN: i32 = (FICT_OFFSET + 15);
/// Release constant `WALDEMATH`.
pub const WALDEMATH: i32 = (FICT_OFFSET + 18);
/// Release constant `PLANET_X_LEVERRIER`.
pub const PLANET_X_LEVERRIER: i32 = (FICT_OFFSET + 11);
/// Release constant `PLANET_X_ADAMS`.
pub const PLANET_X_ADAMS: i32 = (FICT_OFFSET + 12);
/// Release constant `PLANET_X_LOWELL`.
pub const PLANET_X_LOWELL: i32 = (FICT_OFFSET + 13);
/// Release constant `PLANET_X_PICKERING`.
pub const PLANET_X_PICKERING: i32 = (FICT_OFFSET + 14);
/// Release constant `WHITE_MOON`.
pub const WHITE_MOON: i32 = (FICT_OFFSET + 16);
/// Release constant `SELENA`.
pub const SELENA: i32 = WHITE_MOON;
/// Release constant `PROSERPINA`.
pub const PROSERPINA: i32 = (FICT_OFFSET + 17);
/// Release constant `ERIS`.
pub const ERIS: i32 = (136199 + AST_OFFSET);
/// Release constant `SEDNA`.
pub const SEDNA: i32 = (90377 + AST_OFFSET);
/// Release constant `HAUMEA`.
pub const HAUMEA: i32 = (136108 + AST_OFFSET);
/// Release constant `MAKEMAKE`.
pub const MAKEMAKE: i32 = (136472 + AST_OFFSET);
/// Release constant `IXION`.
pub const IXION: i32 = (28978 + AST_OFFSET);
/// Release constant `ORCUS`.
pub const ORCUS: i32 = (90482 + AST_OFFSET);
/// Release constant `QUAOAR`.
pub const QUAOAR: i32 = (50000 + AST_OFFSET);
/// Release constant `NESSUS`.
pub const NESSUS: i32 = (7066 + AST_OFFSET);
/// Release constant `ASBOLUS`.
pub const ASBOLUS: i32 = (8405 + AST_OFFSET);
/// Release constant `CHARIKLO`.
pub const CHARIKLO: i32 = (10199 + AST_OFFSET);
/// Release constant `GONGGONG`.
pub const GONGGONG: i32 = (225088 + AST_OFFSET);
/// Release constant `APOPHIS`.
pub const APOPHIS: i32 = (99942 + AST_OFFSET);
/// Release constant `HYGIEA`.
pub const HYGIEA: i32 = (10 + AST_OFFSET);
/// Release constant `INTERAMNIA`.
pub const INTERAMNIA: i32 = (704 + AST_OFFSET);
/// Release constant `DAVIDA`.
pub const DAVIDA: i32 = (511 + AST_OFFSET);
/// Release constant `EUROPA_AST`.
pub const EUROPA_AST: i32 = (52 + AST_OFFSET);
/// Release constant `SYLVIA`.
pub const SYLVIA: i32 = (87 + AST_OFFSET);
/// Release constant `PSYCHE`.
pub const PSYCHE: i32 = (16 + AST_OFFSET);
/// Release constant `EROS`.
pub const EROS: i32 = (433 + AST_OFFSET);
/// Release constant `AMOR`.
pub const AMOR: i32 = (1221 + AST_OFFSET);
/// Release constant `ICARUS`.
pub const ICARUS: i32 = (1566 + AST_OFFSET);
/// Release constant `TORO`.
pub const TORO: i32 = (1685 + AST_OFFSET);
/// Release constant `SAPPHO`.
pub const SAPPHO: i32 = (80 + AST_OFFSET);
/// Release constant `PANDORA_AST`.
pub const PANDORA_AST: i32 = (55 + AST_OFFSET);
/// Release constant `LILITH_AST`.
pub const LILITH_AST: i32 = (1181 + AST_OFFSET);
/// Release constant `HIDALGO`.
pub const HIDALGO: i32 = (944 + AST_OFFSET);
/// Release constant `TOUTATIS`.
pub const TOUTATIS: i32 = (4179 + AST_OFFSET);
/// Release constant `ITOKAWA`.
pub const ITOKAWA: i32 = (25143 + AST_OFFSET);
/// Release constant `BENNU`.
pub const BENNU: i32 = (101955 + AST_OFFSET);
/// Release constant `RYUGU`.
pub const RYUGU: i32 = (162173 + AST_OFFSET);
/// Release constant `NAIF_ASTEROID_OFFSET`.
pub const NAIF_ASTEROID_OFFSET: i32 = 2000000;
/// Release constant `NAIF_CERES`.
pub const NAIF_CERES: i32 = 2000001;
/// Release constant `NAIF_PALLAS`.
pub const NAIF_PALLAS: i32 = 2000002;
/// Release constant `NAIF_JUNO`.
pub const NAIF_JUNO: i32 = 2000003;
/// Release constant `NAIF_VESTA`.
pub const NAIF_VESTA: i32 = 2000004;
/// Release constant `NAIF_CHIRON`.
pub const NAIF_CHIRON: i32 = 2002060;
/// Release constant `NAIF_PHOLUS`.
pub const NAIF_PHOLUS: i32 = 2005145;
/// Release constant `NAIF_ERIS`.
pub const NAIF_ERIS: i32 = 2136199;
/// Release constant `NAIF_SEDNA`.
pub const NAIF_SEDNA: i32 = 2090377;
/// Release constant `NAIF_HAUMEA`.
pub const NAIF_HAUMEA: i32 = 2136108;
/// Release constant `NAIF_MAKEMAKE`.
pub const NAIF_MAKEMAKE: i32 = 2136472;
/// Release constant `NAIF_IXION`.
pub const NAIF_IXION: i32 = 2028978;
/// Release constant `NAIF_ORCUS`.
pub const NAIF_ORCUS: i32 = 2090482;
/// Release constant `NAIF_QUAOAR`.
pub const NAIF_QUAOAR: i32 = 2050000;
/// Release constant `NAIF_NESSUS`.
pub const NAIF_NESSUS: i32 = 2007066;
/// Release constant `NAIF_ASBOLUS`.
pub const NAIF_ASBOLUS: i32 = 2008405;
/// Release constant `NAIF_CHARIKLO`.
pub const NAIF_CHARIKLO: i32 = 2010199;
/// Release constant `NAIF_GONGGONG`.
pub const NAIF_GONGGONG: i32 = 2225088;
/// Release constant `NAIF_APOPHIS`.
pub const NAIF_APOPHIS: i32 = 2099942;
/// Release constant `NAIF_HYGIEA`.
pub const NAIF_HYGIEA: i32 = 2000010;
/// Release constant `NAIF_EROS`.
pub const NAIF_EROS: i32 = 2000433;
/// Release constant `NAIF_PSYCHE`.
pub const NAIF_PSYCHE: i32 = 2000016;
/// Release constant `NAIF_EUROPA_AST`.
pub const NAIF_EUROPA_AST: i32 = 2000052;
/// Release constant `NAIF_PANDORA_AST`.
pub const NAIF_PANDORA_AST: i32 = 2000055;
/// Release constant `NAIF_SAPPHO`.
pub const NAIF_SAPPHO: i32 = 2000080;
/// Release constant `NAIF_SYLVIA`.
pub const NAIF_SYLVIA: i32 = 2000087;
/// Release constant `NAIF_DAVIDA`.
pub const NAIF_DAVIDA: i32 = 2000511;
/// Release constant `NAIF_INTERAMNIA`.
pub const NAIF_INTERAMNIA: i32 = 2000704;
/// Release constant `NAIF_HIDALGO`.
pub const NAIF_HIDALGO: i32 = 2000944;
/// Release constant `NAIF_LILITH_AST`.
pub const NAIF_LILITH_AST: i32 = 2001181;
/// Release constant `NAIF_AMOR`.
pub const NAIF_AMOR: i32 = 2001221;
/// Release constant `NAIF_ICARUS`.
pub const NAIF_ICARUS: i32 = 2001566;
/// Release constant `NAIF_TORO`.
pub const NAIF_TORO: i32 = 2001685;
/// Release constant `NAIF_TOUTATIS`.
pub const NAIF_TOUTATIS: i32 = 2004179;
/// Release constant `NAIF_ITOKAWA`.
pub const NAIF_ITOKAWA: i32 = 2025143;
/// Release constant `NAIF_BENNU`.
pub const NAIF_BENNU: i32 = 2101955;
/// Release constant `NAIF_RYUGU`.
pub const NAIF_RYUGU: i32 = 2162173;
/// Release constant `FIXSTAR_OFFSET`.
pub const FIXSTAR_OFFSET: i32 = 1000000;
/// Release constant `REGULUS`.
pub const REGULUS: i32 = (FIXSTAR_OFFSET + 1);
/// Release constant `SPICA_STAR`.
pub const SPICA_STAR: i32 = (FIXSTAR_OFFSET + 2);
/// Release constant `ALGOL`.
pub const ALGOL: i32 = (FIXSTAR_OFFSET + 3);
/// Release constant `SIRIUS`.
pub const SIRIUS: i32 = (FIXSTAR_OFFSET + 4);
/// Release constant `ALDEBARAN`.
pub const ALDEBARAN: i32 = (FIXSTAR_OFFSET + 5);
/// Release constant `ANTARES`.
pub const ANTARES: i32 = (FIXSTAR_OFFSET + 6);
/// Release constant `VEGA`.
pub const VEGA: i32 = (FIXSTAR_OFFSET + 7);
/// Release constant `POLARIS`.
pub const POLARIS: i32 = (FIXSTAR_OFFSET + 8);
/// Release constant `FOMALHAUT`.
pub const FOMALHAUT: i32 = (FIXSTAR_OFFSET + 9);
/// Release constant `BETELGEUSE`.
pub const BETELGEUSE: i32 = (FIXSTAR_OFFSET + 10);
/// Release constant `RIGEL`.
pub const RIGEL: i32 = (FIXSTAR_OFFSET + 11);
/// Release constant `PROCYON`.
pub const PROCYON: i32 = (FIXSTAR_OFFSET + 12);
/// Release constant `CAPELLA`.
pub const CAPELLA: i32 = (FIXSTAR_OFFSET + 13);
/// Release constant `ARCTURUS`.
pub const ARCTURUS: i32 = (FIXSTAR_OFFSET + 14);
/// Release constant `DENEB`.
pub const DENEB: i32 = (FIXSTAR_OFFSET + 15);
/// Release constant `POLLUX`.
pub const POLLUX: i32 = (FIXSTAR_OFFSET + 16);
/// Release constant `CASTOR`.
pub const CASTOR: i32 = (FIXSTAR_OFFSET + 17);
/// Release constant `ALTAIR`.
pub const ALTAIR: i32 = (FIXSTAR_OFFSET + 18);
/// Release constant `ACHERNAR`.
pub const ACHERNAR: i32 = (FIXSTAR_OFFSET + 19);
/// Release constant `CANOPUS`.
pub const CANOPUS: i32 = (FIXSTAR_OFFSET + 20);
/// Release constant `ACRUX`.
pub const ACRUX: i32 = (FIXSTAR_OFFSET + 21);
/// Release constant `MIMOSA`.
pub const MIMOSA: i32 = (FIXSTAR_OFFSET + 22);
/// Release constant `GACRUX`.
pub const GACRUX: i32 = (FIXSTAR_OFFSET + 23);
/// Release constant `HADAR`.
pub const HADAR: i32 = (FIXSTAR_OFFSET + 24);
/// Release constant `RIGIL_KENT`.
pub const RIGIL_KENT: i32 = (FIXSTAR_OFFSET + 25);
/// Release constant `SHAULA`.
pub const SHAULA: i32 = (FIXSTAR_OFFSET + 26);
/// Release constant `BELLATRIX`.
pub const BELLATRIX: i32 = (FIXSTAR_OFFSET + 27);
/// Release constant `ELNATH`.
pub const ELNATH: i32 = (FIXSTAR_OFFSET + 28);
/// Release constant `MIRA`.
pub const MIRA: i32 = (FIXSTAR_OFFSET + 29);
/// Release constant `ALNILAM`.
pub const ALNILAM: i32 = (FIXSTAR_OFFSET + 30);
/// Release constant `ALNITAK`.
pub const ALNITAK: i32 = (FIXSTAR_OFFSET + 31);
/// Release constant `MINTAKA`.
pub const MINTAKA: i32 = (FIXSTAR_OFFSET + 32);
/// Release constant `SAIPH`.
pub const SAIPH: i32 = (FIXSTAR_OFFSET + 33);
/// Release constant `DIPHDA`.
pub const DIPHDA: i32 = (FIXSTAR_OFFSET + 34);
/// Release constant `ALPHARD`.
pub const ALPHARD: i32 = (FIXSTAR_OFFSET + 35);
/// Release constant `RASALHAGUE`.
pub const RASALHAGUE: i32 = (FIXSTAR_OFFSET + 36);
/// Release constant `ETAMIN`.
pub const ETAMIN: i32 = (FIXSTAR_OFFSET + 37);
/// Release constant `KOCHAB`.
pub const KOCHAB: i32 = (FIXSTAR_OFFSET + 38);
/// Release constant `ALKAID`.
pub const ALKAID: i32 = (FIXSTAR_OFFSET + 39);
/// Release constant `DUBHE`.
pub const DUBHE: i32 = (FIXSTAR_OFFSET + 40);
/// Release constant `MERAK`.
pub const MERAK: i32 = (FIXSTAR_OFFSET + 41);
/// Release constant `ALIOTH`.
pub const ALIOTH: i32 = (FIXSTAR_OFFSET + 42);
/// Release constant `MIZAR`.
pub const MIZAR: i32 = (FIXSTAR_OFFSET + 43);
/// Release constant `ALCOR`.
pub const ALCOR: i32 = (FIXSTAR_OFFSET + 44);
/// Release constant `VINDEMIATRIX`.
pub const VINDEMIATRIX: i32 = (FIXSTAR_OFFSET + 45);
/// Release constant `ZUBENELGENUBI`.
pub const ZUBENELGENUBI: i32 = (FIXSTAR_OFFSET + 46);
/// Release constant `ZUBENESCHAMALI`.
pub const ZUBENESCHAMALI: i32 = (FIXSTAR_OFFSET + 47);
/// Release constant `UNUKALHAI`.
pub const UNUKALHAI: i32 = (FIXSTAR_OFFSET + 48);
/// Release constant `ALGIEBA`.
pub const ALGIEBA: i32 = (FIXSTAR_OFFSET + 49);
/// Release constant `DENEBOLA`.
pub const DENEBOLA: i32 = (FIXSTAR_OFFSET + 50);
/// Release constant `MARKAB`.
pub const MARKAB: i32 = (FIXSTAR_OFFSET + 51);
/// Release constant `SCHEAT`.
pub const SCHEAT: i32 = (FIXSTAR_OFFSET + 52);
/// Release constant `ALCYONE`.
pub const ALCYONE: i32 = (FIXSTAR_OFFSET + 53);
/// Release constant `ALGORAB`.
pub const ALGORAB: i32 = (FIXSTAR_OFFSET + 54);
/// Release constant `ALPHECCA`.
pub const ALPHECCA: i32 = (FIXSTAR_OFFSET + 55);
/// Release constant `DENEB_ALGEDI`.
pub const DENEB_ALGEDI: i32 = (FIXSTAR_OFFSET + 56);
/// Release constant `ASTEROPE`.
pub const ASTEROPE: i32 = (FIXSTAR_OFFSET + 57);
/// Release constant `CELAENO`.
pub const CELAENO: i32 = (FIXSTAR_OFFSET + 58);
/// Release constant `ELECTRA`.
pub const ELECTRA: i32 = (FIXSTAR_OFFSET + 59);
/// Release constant `MAIA`.
pub const MAIA: i32 = (FIXSTAR_OFFSET + 60);
/// Release constant `MEROPE`.
pub const MEROPE: i32 = (FIXSTAR_OFFSET + 61);
/// Release constant `TAYGETA`.
pub const TAYGETA: i32 = (FIXSTAR_OFFSET + 62);
/// Release constant `ATLAS`.
pub const ATLAS: i32 = (FIXSTAR_OFFSET + 63);
/// Release constant `PLEIONE`.
pub const PLEIONE: i32 = (FIXSTAR_OFFSET + 64);
/// Release constant `PRIMA_HYADUM`.
pub const PRIMA_HYADUM: i32 = (FIXSTAR_OFFSET + 65);
/// Release constant `SECUNDA_HYADUM`.
pub const SECUNDA_HYADUM: i32 = (FIXSTAR_OFFSET + 66);
/// Release constant `THETA_TAURI`.
pub const THETA_TAURI: i32 = (FIXSTAR_OFFSET + 67);
/// Release constant `AIN`.
pub const AIN: i32 = (FIXSTAR_OFFSET + 68);
/// Release constant `MEISSA`.
pub const MEISSA: i32 = (FIXSTAR_OFFSET + 69);
/// Release constant `PHECDA`.
pub const PHECDA: i32 = (FIXSTAR_OFFSET + 70);
/// Release constant `MEGREZ`.
pub const MEGREZ: i32 = (FIXSTAR_OFFSET + 71);
/// Release constant `DELTA_CRUCIS`.
pub const DELTA_CRUCIS: i32 = (FIXSTAR_OFFSET + 72);
/// Release constant `MENKENT`.
pub const MENKENT: i32 = (FIXSTAR_OFFSET + 73);
/// Release constant `MUHLIFAIN`.
pub const MUHLIFAIN: i32 = (FIXSTAR_OFFSET + 74);
/// Release constant `EPSILON_CENTAURI`.
pub const EPSILON_CENTAURI: i32 = (FIXSTAR_OFFSET + 75);
/// Release constant `ETA_CENTAURI`.
pub const ETA_CENTAURI: i32 = (FIXSTAR_OFFSET + 76);
/// Release constant `ZETA_CENTAURI`.
pub const ZETA_CENTAURI: i32 = (FIXSTAR_OFFSET + 77);
/// Release constant `SARGAS`.
pub const SARGAS: i32 = (FIXSTAR_OFFSET + 78);
/// Release constant `DSCHUBBA`.
pub const DSCHUBBA: i32 = (FIXSTAR_OFFSET + 79);
/// Release constant `GRAFFIAS`.
pub const GRAFFIAS: i32 = (FIXSTAR_OFFSET + 80);
/// Release constant `LESATH`.
pub const LESATH: i32 = (FIXSTAR_OFFSET + 81);
/// Release constant `ZOSMA`.
pub const ZOSMA: i32 = (FIXSTAR_OFFSET + 82);
/// Release constant `HAMAL`.
pub const HAMAL: i32 = (FIXSTAR_OFFSET + 83);
/// Release constant `SHERATAN`.
pub const SHERATAN: i32 = (FIXSTAR_OFFSET + 84);
/// Release constant `MESARTHIM`.
pub const MESARTHIM: i32 = (FIXSTAR_OFFSET + 85);
/// Release constant `ACUBENS`.
pub const ACUBENS: i32 = (FIXSTAR_OFFSET + 86);
/// Release constant `TARF`.
pub const TARF: i32 = (FIXSTAR_OFFSET + 87);
/// Release constant `ASELLUS_BOREALIS`.
pub const ASELLUS_BOREALIS: i32 = (FIXSTAR_OFFSET + 88);
/// Release constant `ASELLUS_AUSTRALIS`.
pub const ASELLUS_AUSTRALIS: i32 = (FIXSTAR_OFFSET + 89);
/// Release constant `KAUS_AUSTRALIS`.
pub const KAUS_AUSTRALIS: i32 = (FIXSTAR_OFFSET + 90);
/// Release constant `NUNKI`.
pub const NUNKI: i32 = (FIXSTAR_OFFSET + 91);
/// Release constant `KAUS_MEDIA`.
pub const KAUS_MEDIA: i32 = (FIXSTAR_OFFSET + 92);
/// Release constant `KAUS_BOREALIS`.
pub const KAUS_BOREALIS: i32 = (FIXSTAR_OFFSET + 93);
/// Release constant `ASCELLA`.
pub const ASCELLA: i32 = (FIXSTAR_OFFSET + 94);
/// Release constant `ALGEDI`.
pub const ALGEDI: i32 = (FIXSTAR_OFFSET + 95);
/// Release constant `DABIH`.
pub const DABIH: i32 = (FIXSTAR_OFFSET + 96);
/// Release constant `NASHIRA`.
pub const NASHIRA: i32 = (FIXSTAR_OFFSET + 97);
/// Release constant `SADALSUUD`.
pub const SADALSUUD: i32 = (FIXSTAR_OFFSET + 98);
/// Release constant `SADALMELIK`.
pub const SADALMELIK: i32 = (FIXSTAR_OFFSET + 99);
/// Release constant `SKAT`.
pub const SKAT: i32 = (FIXSTAR_OFFSET + 100);
/// Release constant `ETA_PISCIUM`.
pub const ETA_PISCIUM: i32 = (FIXSTAR_OFFSET + 101);
/// Release constant `ALRESCHA`.
pub const ALRESCHA: i32 = (FIXSTAR_OFFSET + 102);
/// Release constant `ALPHERATZ`.
pub const ALPHERATZ: i32 = (FIXSTAR_OFFSET + 103);
/// Release constant `ALGENIB`.
pub const ALGENIB: i32 = (FIXSTAR_OFFSET + 104);
/// Release constant `PROPUS`.
pub const PROPUS: i32 = (FIXSTAR_OFFSET + 105);
/// Release constant `ALHENA`.
pub const ALHENA: i32 = (FIXSTAR_OFFSET + 106);
/// Release constant `TEJAT`.
pub const TEJAT: i32 = (FIXSTAR_OFFSET + 107);
/// Release constant `WASAT`.
pub const WASAT: i32 = (FIXSTAR_OFFSET + 108);
/// Release constant `ADHARA`.
pub const ADHARA: i32 = (FIXSTAR_OFFSET + 109);
/// Release constant `WEZEN`.
pub const WEZEN: i32 = (FIXSTAR_OFFSET + 110);
/// Release constant `THUBAN`.
pub const THUBAN: i32 = (FIXSTAR_OFFSET + 111);
/// Release constant `RASALGETHI`.
pub const RASALGETHI: i32 = (FIXSTAR_OFFSET + 112);
/// Release constant `ALBIREO`.
pub const ALBIREO: i32 = (FIXSTAR_OFFSET + 113);
/// Release constant `MIRACH`.
pub const MIRACH: i32 = (FIXSTAR_OFFSET + 114);
/// Release constant `ALMACH`.
pub const ALMACH: i32 = (FIXSTAR_OFFSET + 115);
/// Release constant `MENKAR`.
pub const MENKAR: i32 = (FIXSTAR_OFFSET + 116);
/// Release constant `ANGLE_OFFSET`.
pub const ANGLE_OFFSET: i32 = 2000000;
/// Release constant `ASCENDANT`.
pub const ASCENDANT: i32 = (ANGLE_OFFSET + 1);
/// Release constant `DESCENDANT`.
pub const DESCENDANT: i32 = (ANGLE_OFFSET + 3);
/// Release constant `IC`.
pub const IC: i32 = (ANGLE_OFFSET + 4);
/// Release constant `ANTIVERTEX`.
pub const ANTIVERTEX: i32 = (ANGLE_OFFSET + 6);
/// Release constant `ARABIC_OFFSET`.
pub const ARABIC_OFFSET: i32 = 3000000;
/// Release constant `PARS_FORTUNAE`.
pub const PARS_FORTUNAE: i32 = (ARABIC_OFFSET + 1);
/// Release constant `PARS_SPIRITUS`.
pub const PARS_SPIRITUS: i32 = (ARABIC_OFFSET + 2);
/// Release constant `PARS_AMORIS`.
pub const PARS_AMORIS: i32 = (ARABIC_OFFSET + 3);
/// Release constant `PARS_FIDEI`.
pub const PARS_FIDEI: i32 = (ARABIC_OFFSET + 4);
/// Release constant `FLG_JPLEPH`.
pub const FLG_JPLEPH: i32 = 1;
/// Release constant `FLG_SWIEPH`.
pub const FLG_SWIEPH: i32 = 2;
/// Release constant `FLG_MOSEPH`.
pub const FLG_MOSEPH: i32 = 4;
/// Release constant `FLG_HELCTR`.
pub const FLG_HELCTR: i32 = 8;
/// Release constant `FLG_TRUEPOS`.
pub const FLG_TRUEPOS: i32 = 16;
/// Release constant `FLG_J2000`.
pub const FLG_J2000: i32 = 32;
/// Release constant `FLG_NONUT`.
pub const FLG_NONUT: i32 = 64;
/// Release constant `FLG_SPEED3`.
pub const FLG_SPEED3: i32 = 128;
/// Release constant `FLG_SPEED`.
pub const FLG_SPEED: i32 = 256;
/// Release constant `FLG_NOGDEFL`.
pub const FLG_NOGDEFL: i32 = 512;
/// Release constant `FLG_NOABERR`.
pub const FLG_NOABERR: i32 = 1024;
/// Release constant `FLG_ASTROMETRIC`.
pub const FLG_ASTROMETRIC: i32 = (FLG_NOABERR | FLG_NOGDEFL);
/// Release constant `FLG_EQUATORIAL`.
pub const FLG_EQUATORIAL: i32 = 2048;
/// Release constant `FLG_XYZ`.
pub const FLG_XYZ: i32 = 4096;
/// Release constant `FLG_RADIANS`.
pub const FLG_RADIANS: i32 = 8192;
/// Release constant `FLG_BARYCTR`.
pub const FLG_BARYCTR: i32 = 16384;
/// Release constant `FLG_TOPOCTR`.
pub const FLG_TOPOCTR: i32 = 32768;
/// Release constant `FLG_SIDEREAL`.
pub const FLG_SIDEREAL: i32 = 65536;
/// Release constant `FLG_ICRS`.
pub const FLG_ICRS: i32 = 131072;
/// Release constant `SIDM_FAGAN_BRADLEY`.
pub const SIDM_FAGAN_BRADLEY: i32 = 0;
/// Release constant `SIDM_LAHIRI`.
pub const SIDM_LAHIRI: i32 = 1;
/// Release constant `SIDM_DELUCE`.
pub const SIDM_DELUCE: i32 = 2;
/// Release constant `SIDM_RAMAN`.
pub const SIDM_RAMAN: i32 = 3;
/// Release constant `SIDM_USHASHASHI`.
pub const SIDM_USHASHASHI: i32 = 4;
/// Release constant `SIDM_KRISHNAMURTI`.
pub const SIDM_KRISHNAMURTI: i32 = 5;
/// Release constant `SIDM_DJWHAL_KHUL`.
pub const SIDM_DJWHAL_KHUL: i32 = 6;
/// Release constant `SIDM_YUKTESHWAR`.
pub const SIDM_YUKTESHWAR: i32 = 7;
/// Release constant `SIDM_JN_BHASIN`.
pub const SIDM_JN_BHASIN: i32 = 8;
/// Release constant `SIDM_BABYL_KUGLER1`.
pub const SIDM_BABYL_KUGLER1: i32 = 9;
/// Release constant `SIDM_BABYL_KUGLER2`.
pub const SIDM_BABYL_KUGLER2: i32 = 10;
/// Release constant `SIDM_BABYL_KUGLER3`.
pub const SIDM_BABYL_KUGLER3: i32 = 11;
/// Release constant `SIDM_BABYL_HUBER`.
pub const SIDM_BABYL_HUBER: i32 = 12;
/// Release constant `SIDM_BABYL_ETPSC`.
pub const SIDM_BABYL_ETPSC: i32 = 13;
/// Release constant `SIDM_BABYL_BRITTON`.
pub const SIDM_BABYL_BRITTON: i32 = 38;
/// Release constant `SIDM_ALDEBARAN_15TAU`.
pub const SIDM_ALDEBARAN_15TAU: i32 = 14;
/// Release constant `SIDM_TRUE_CITRA`.
pub const SIDM_TRUE_CITRA: i32 = 27;
/// Release constant `SIDM_TRUE_REVATI`.
pub const SIDM_TRUE_REVATI: i32 = 28;
/// Release constant `SIDM_TRUE_PUSHYA`.
pub const SIDM_TRUE_PUSHYA: i32 = 29;
/// Release constant `SIDM_TRUE_MULA`.
pub const SIDM_TRUE_MULA: i32 = 35;
/// Release constant `SIDM_TRUE_SHEORAN`.
pub const SIDM_TRUE_SHEORAN: i32 = 39;
/// Release constant `SIDM_HIPPARCHOS`.
pub const SIDM_HIPPARCHOS: i32 = 15;
/// Release constant `SIDM_SASSANIAN`.
pub const SIDM_SASSANIAN: i32 = 16;
/// Release constant `SIDM_J2000`.
pub const SIDM_J2000: i32 = 18;
/// Release constant `SIDM_J1900`.
pub const SIDM_J1900: i32 = 19;
/// Release constant `SIDM_B1950`.
pub const SIDM_B1950: i32 = 20;
/// Release constant `SIDM_SURYASIDDHANTA`.
pub const SIDM_SURYASIDDHANTA: i32 = 21;
/// Release constant `SIDM_SURYASIDDHANTA_MSUN`.
pub const SIDM_SURYASIDDHANTA_MSUN: i32 = 22;
/// Release constant `SIDM_ARYABHATA`.
pub const SIDM_ARYABHATA: i32 = 23;
/// Release constant `SIDM_ARYABHATA_MSUN`.
pub const SIDM_ARYABHATA_MSUN: i32 = 24;
/// Release constant `SIDM_ARYABHATA_522`.
pub const SIDM_ARYABHATA_522: i32 = 37;
/// Release constant `SIDM_SS_REVATI`.
pub const SIDM_SS_REVATI: i32 = 25;
/// Release constant `SIDM_SS_CITRA`.
pub const SIDM_SS_CITRA: i32 = 26;
/// Release constant `SIDM_GALCENT_0SAG`.
pub const SIDM_GALCENT_0SAG: i32 = 17;
/// Release constant `SIDM_GALCENT_RGILBRAND`.
pub const SIDM_GALCENT_RGILBRAND: i32 = 30;
/// Release constant `SIDM_GALCENT_MULA_WILHELM`.
pub const SIDM_GALCENT_MULA_WILHELM: i32 = 36;
/// Release constant `SIDM_GALCENT_COCHRANE`.
pub const SIDM_GALCENT_COCHRANE: i32 = 40;
/// Release constant `SIDM_GALEQU_IAU1958`.
pub const SIDM_GALEQU_IAU1958: i32 = 31;
/// Release constant `SIDM_GALEQU_TRUE`.
pub const SIDM_GALEQU_TRUE: i32 = 32;
/// Release constant `SIDM_GALEQU_MULA`.
pub const SIDM_GALEQU_MULA: i32 = 33;
/// Release constant `SIDM_GALEQU_FIORENZA`.
pub const SIDM_GALEQU_FIORENZA: i32 = 41;
/// Release constant `SIDM_GALALIGN_MARDYKS`.
pub const SIDM_GALALIGN_MARDYKS: i32 = 34;
/// Release constant `SIDM_VALENS_MOON`.
pub const SIDM_VALENS_MOON: i32 = 42;
/// Release constant `SIDM_LAHIRI_1940`.
pub const SIDM_LAHIRI_1940: i32 = 43;
/// Release constant `SIDM_LAHIRI_VP285`.
pub const SIDM_LAHIRI_VP285: i32 = 44;
/// Release constant `SIDM_KRISHNAMURTI_VP291`.
pub const SIDM_KRISHNAMURTI_VP291: i32 = 45;
/// Release constant `SIDM_LAHIRI_ICRC`.
pub const SIDM_LAHIRI_ICRC: i32 = 46;
/// Release constant `SIDM_USER`.
pub const SIDM_USER: i32 = 255;
/// Release constant `JUL_CAL`.
pub const JUL_CAL: i32 = 0;
/// Release constant `GREG_CAL`.
pub const GREG_CAL: i32 = 1;
/// Release constant `AUNIT` (bit-exact).
pub const AUNIT: f64 = f64::from_bits(0x41A1D55D5D666666);
/// Release constant `ECL_CENTRAL`.
pub const ECL_CENTRAL: i32 = 1;
/// Release constant `ECL_NONCENTRAL`.
pub const ECL_NONCENTRAL: i32 = 2;
/// Release constant `ECL_TOTAL`.
pub const ECL_TOTAL: i32 = 4;
/// Release constant `ECL_ANNULAR`.
pub const ECL_ANNULAR: i32 = 8;
/// Release constant `ECL_PARTIAL`.
pub const ECL_PARTIAL: i32 = 16;
/// Release constant `ECL_ANNULAR_TOTAL`.
pub const ECL_ANNULAR_TOTAL: i32 = 32;
/// Release constant `ECL_PENUMBRAL`.
pub const ECL_PENUMBRAL: i32 = 64;
/// Release constant `ECL_GRAZING`.
pub const ECL_GRAZING: i32 = 65536;
/// Release constant `ECL_ALLTYPES_SOLAR`.
pub const ECL_ALLTYPES_SOLAR: i32 =
    (((((ECL_CENTRAL | ECL_NONCENTRAL) | ECL_TOTAL) | ECL_ANNULAR) | ECL_PARTIAL)
        | ECL_ANNULAR_TOTAL);
/// Release constant `ECL_ALLTYPES_LUNAR`.
pub const ECL_ALLTYPES_LUNAR: i32 = ((ECL_TOTAL | ECL_PARTIAL) | ECL_PENUMBRAL);
/// Release constant `ECL_VISIBLE`.
pub const ECL_VISIBLE: i32 = 128;
/// Release constant `ECL_MAX_VISIBLE`.
pub const ECL_MAX_VISIBLE: i32 = 256;
/// Release constant `ECL_1ST_VISIBLE`.
pub const ECL_1ST_VISIBLE: i32 = 512;
/// Release constant `ECL_2ND_VISIBLE`.
pub const ECL_2ND_VISIBLE: i32 = 1024;
/// Release constant `ECL_3RD_VISIBLE`.
pub const ECL_3RD_VISIBLE: i32 = 2048;
/// Release constant `ECL_4TH_VISIBLE`.
pub const ECL_4TH_VISIBLE: i32 = 4096;
/// Release constant `ECL_ONE_TRY`.
pub const ECL_ONE_TRY: i32 = 32768;
/// Release constant `NODBIT_MEAN`.
pub const NODBIT_MEAN: i32 = 1;
/// Release constant `NODBIT_OSCU`.
pub const NODBIT_OSCU: i32 = 2;
/// Release constant `NODBIT_OSCU_BAR`.
pub const NODBIT_OSCU_BAR: i32 = 4;
/// Release constant `NODBIT_FOPOINT`.
pub const NODBIT_FOPOINT: i32 = 256;
/// Release constant `CALC_RISE`.
pub const CALC_RISE: i32 = 1;
/// Release constant `CALC_SET`.
pub const CALC_SET: i32 = 2;
/// Release constant `CALC_MTRANSIT`.
pub const CALC_MTRANSIT: i32 = 4;
/// Release constant `CALC_ITRANSIT`.
pub const CALC_ITRANSIT: i32 = 8;
/// Release constant `BIT_DISC_CENTER`.
pub const BIT_DISC_CENTER: i32 = 256;
/// Release constant `BIT_DISC_BOTTOM`.
pub const BIT_DISC_BOTTOM: i32 = 8192;
/// Release constant `BIT_NO_REFRACTION`.
pub const BIT_NO_REFRACTION: i32 = 512;
/// Release constant `BIT_CIVIL_TWILIGHT`.
pub const BIT_CIVIL_TWILIGHT: i32 = 1024;
/// Release constant `BIT_NAUTIC_TWILIGHT`.
pub const BIT_NAUTIC_TWILIGHT: i32 = 2048;
/// Release constant `BIT_ASTRO_TWILIGHT`.
pub const BIT_ASTRO_TWILIGHT: i32 = 4096;
/// Release constant `BIT_FIXED_DISC_SIZE`.
pub const BIT_FIXED_DISC_SIZE: i32 = 16384;
/// Release constant `HELIACAL_RISING`.
pub const HELIACAL_RISING: i32 = 1;
/// Release constant `HELIACAL_SETTING`.
pub const HELIACAL_SETTING: i32 = 2;
/// Release constant `MORNING_FIRST`.
pub const MORNING_FIRST: i32 = HELIACAL_RISING;
/// Release constant `EVENING_LAST`.
pub const EVENING_LAST: i32 = HELIACAL_SETTING;
/// Release constant `EVENING_FIRST`.
pub const EVENING_FIRST: i32 = 3;
/// Release constant `MORNING_LAST`.
pub const MORNING_LAST: i32 = 4;
/// Release constant `HELFLAG_LONG_SEARCH`.
pub const HELFLAG_LONG_SEARCH: i32 = 128;
/// Release constant `HELFLAG_HIGH_PRECISION`.
pub const HELFLAG_HIGH_PRECISION: i32 = 256;
/// Release constant `HELFLAG_OPTICAL_PARAMS`.
pub const HELFLAG_OPTICAL_PARAMS: i32 = (1 << 9);
/// Release constant `HELFLAG_NO_DETAILS`.
pub const HELFLAG_NO_DETAILS: i32 = (1 << 10);
/// Release constant `HELFLAG_SEARCH_1_PERIOD`.
pub const HELFLAG_SEARCH_1_PERIOD: i32 = (1 << 11);
/// Release constant `HELFLAG_VISLIM_DARK`.
pub const HELFLAG_VISLIM_DARK: i32 = (1 << 12);
/// Release constant `HELFLAG_VISLIM_NOMOON`.
pub const HELFLAG_VISLIM_NOMOON: i32 = (1 << 13);
/// Release constant `HELFLAG_VISLIM_PHOTOPIC`.
pub const HELFLAG_VISLIM_PHOTOPIC: i32 = (1 << 14);
/// Release constant `HELFLAG_VISLIM_SCOTOPIC`.
pub const HELFLAG_VISLIM_SCOTOPIC: i32 = (1 << 15);
/// Release constant `HELFLAG_AV`.
pub const HELFLAG_AV: i32 = (1 << 16);
/// Release constant `HELFLAG_AVKIND_VR`.
pub const HELFLAG_AVKIND_VR: i32 = (1 << 16);
/// Release constant `HELFLAG_AVKIND_PTO`.
pub const HELFLAG_AVKIND_PTO: i32 = (1 << 17);
/// Release constant `HELFLAG_AVKIND_MIN7`.
pub const HELFLAG_AVKIND_MIN7: i32 = (1 << 18);
/// Release constant `HELFLAG_AVKIND_MIN9`.
pub const HELFLAG_AVKIND_MIN9: i32 = (1 << 19);
/// Release constant `HELFLAG_AVKIND`.
pub const HELFLAG_AVKIND: i32 =
    (((HELFLAG_AVKIND_VR | HELFLAG_AVKIND_PTO) | HELFLAG_AVKIND_MIN7) | HELFLAG_AVKIND_MIN9);
/// Release constant `HELFLAG_BELOW_HORIZON`.
pub const HELFLAG_BELOW_HORIZON: i32 = -2;
/// Release constant `HELFLAG_PHOTOPIC`.
pub const HELFLAG_PHOTOPIC: i32 = 0;
/// Release constant `HELFLAG_SCOTOPIC`.
pub const HELFLAG_SCOTOPIC: i32 = 1;
/// Release constant `HELFLAG_MIXED`.
pub const HELFLAG_MIXED: i32 = 2;
/// Release constant `SPLIT_DEG_ROUND_SEC`.
pub const SPLIT_DEG_ROUND_SEC: i32 = 1;
/// Release constant `SPLIT_DEG_ROUND_MIN`.
pub const SPLIT_DEG_ROUND_MIN: i32 = 2;
/// Release constant `SPLIT_DEG_ROUND_DEG`.
pub const SPLIT_DEG_ROUND_DEG: i32 = 4;
/// Release constant `SPLIT_DEG_ZODIACAL`.
pub const SPLIT_DEG_ZODIACAL: i32 = 8;
/// Release constant `SPLIT_DEG_NAKSHATRA`.
pub const SPLIT_DEG_NAKSHATRA: i32 = 1024;
/// Release constant `SPLIT_DEG_KEEP_SIGN`.
pub const SPLIT_DEG_KEEP_SIGN: i32 = 16;
/// Release constant `SPLIT_DEG_KEEP_DEG`.
pub const SPLIT_DEG_KEEP_DEG: i32 = 32;
/// Release constant `TIDAL_DE200` (bit-exact).
pub const TIDAL_DE200: f64 = f64::from_bits(0xC037E504816F0069);
/// Release constant `TIDAL_DE403` (bit-exact).
pub const TIDAL_DE403: f64 = f64::from_bits(0xC039947AE147AE14);
/// Release constant `TIDAL_DE404` (bit-exact).
pub const TIDAL_DE404: f64 = f64::from_bits(0xC039947AE147AE14);
/// Release constant `TIDAL_DE405` (bit-exact).
pub const TIDAL_DE405: f64 = f64::from_bits(0xC039D374BC6A7EFA);
/// Release constant `TIDAL_DE406` (bit-exact).
pub const TIDAL_DE406: f64 = f64::from_bits(0xC039D374BC6A7EFA);
/// Release constant `TIDAL_DE421` (bit-exact).
pub const TIDAL_DE421: f64 = f64::from_bits(0xC039D9999999999A);
/// Release constant `TIDAL_DE422` (bit-exact).
pub const TIDAL_DE422: f64 = f64::from_bits(0xC039D9999999999A);
/// Release constant `TIDAL_DE430` (bit-exact).
pub const TIDAL_DE430: f64 = f64::from_bits(0xC039D1EB851EB852);
/// Release constant `TIDAL_DE431` (bit-exact).
pub const TIDAL_DE431: f64 = f64::from_bits(0xC039CCCCCCCCCCCD);
/// Release constant `TIDAL_DE440` (bit-exact).
pub const TIDAL_DE440: f64 = f64::from_bits(0xC039EF9DB22D0E56);
/// Release constant `TIDAL_DE441` (bit-exact).
pub const TIDAL_DE441: f64 = f64::from_bits(0xC039EF9DB22D0E56);
/// Release constant `TIDAL_DEFAULT` (bit-exact).
pub const TIDAL_DEFAULT: f64 = f64::from_bits(0xC039CCCCCCCCCCCD);
/// Release constant `TIDAL_AUTOMATIC`.
pub const TIDAL_AUTOMATIC: i32 = 999999;
/// Release constant `J2000` (bit-exact).
pub const J2000: f64 = f64::from_bits(0x4142B42C80000000);
/// Release constant `J1991_25` (bit-exact).
pub const J1991_25: f64 = f64::from_bits(0x4142ADEE88000000);
/// Release constant `J1900` (bit-exact).
pub const J1900: f64 = f64::from_bits(0x41426CD600000000);
/// Release constant `B1950` (bit-exact).
pub const B1950: f64 = f64::from_bits(0x414290813633E7F9);
/// Release constant `DAYS_PER_JULIAN_YEAR` (bit-exact).
pub const DAYS_PER_JULIAN_YEAR: f64 = f64::from_bits(0x4076D40000000000);
/// Release constant `DAYS_PER_JULIAN_CENTURY` (bit-exact).
pub const DAYS_PER_JULIAN_CENTURY: f64 = f64::from_bits(0x40E1D5A000000000);
/// Release constant `MOON_OFFSET`.
pub const MOON_OFFSET: i32 = 9000;
/// Release constant `MOON_IO`.
pub const MOON_IO: i32 = (MOON_OFFSET + 1);
/// Release constant `MOON_EUROPA`.
pub const MOON_EUROPA: i32 = (MOON_OFFSET + 2);
/// Release constant `MOON_GANYMEDE`.
pub const MOON_GANYMEDE: i32 = (MOON_OFFSET + 3);
/// Release constant `MOON_CALLISTO`.
pub const MOON_CALLISTO: i32 = (MOON_OFFSET + 4);
/// Release constant `MOON_MIMAS`.
pub const MOON_MIMAS: i32 = (MOON_OFFSET + 11);
/// Release constant `MOON_ENCELADUS`.
pub const MOON_ENCELADUS: i32 = (MOON_OFFSET + 12);
/// Release constant `MOON_TETHYS`.
pub const MOON_TETHYS: i32 = (MOON_OFFSET + 13);
/// Release constant `MOON_DIONE`.
pub const MOON_DIONE: i32 = (MOON_OFFSET + 14);
/// Release constant `MOON_RHEA`.
pub const MOON_RHEA: i32 = (MOON_OFFSET + 15);
/// Release constant `MOON_TITAN`.
pub const MOON_TITAN: i32 = (MOON_OFFSET + 16);
/// Release constant `MOON_HYPERION`.
pub const MOON_HYPERION: i32 = (MOON_OFFSET + 17);
/// Release constant `MOON_IAPETUS`.
pub const MOON_IAPETUS: i32 = (MOON_OFFSET + 18);
/// Release constant `MOON_MIRANDA`.
pub const MOON_MIRANDA: i32 = (MOON_OFFSET + 21);
/// Release constant `MOON_ARIEL`.
pub const MOON_ARIEL: i32 = (MOON_OFFSET + 22);
/// Release constant `MOON_UMBRIEL`.
pub const MOON_UMBRIEL: i32 = (MOON_OFFSET + 23);
/// Release constant `MOON_TITANIA`.
pub const MOON_TITANIA: i32 = (MOON_OFFSET + 24);
/// Release constant `MOON_OBERON`.
pub const MOON_OBERON: i32 = (MOON_OFFSET + 25);
/// Release constant `MOON_TRITON`.
pub const MOON_TRITON: i32 = (MOON_OFFSET + 31);
/// Release constant `MOON_PHOBOS`.
pub const MOON_PHOBOS: i32 = (MOON_OFFSET + 41);
/// Release constant `MOON_DEIMOS`.
pub const MOON_DEIMOS: i32 = (MOON_OFFSET + 42);
/// Release constant `MOON_CHARON`.
pub const MOON_CHARON: i32 = (MOON_OFFSET + 51);
/// Release constant `NAIF_IO`.
pub const NAIF_IO: i32 = 501;
/// Release constant `NAIF_EUROPA`.
pub const NAIF_EUROPA: i32 = 502;
/// Release constant `NAIF_GANYMEDE`.
pub const NAIF_GANYMEDE: i32 = 503;
/// Release constant `NAIF_CALLISTO`.
pub const NAIF_CALLISTO: i32 = 504;
/// Release constant `NAIF_MIMAS`.
pub const NAIF_MIMAS: i32 = 601;
/// Release constant `NAIF_ENCELADUS`.
pub const NAIF_ENCELADUS: i32 = 602;
/// Release constant `NAIF_TETHYS`.
pub const NAIF_TETHYS: i32 = 603;
/// Release constant `NAIF_DIONE`.
pub const NAIF_DIONE: i32 = 604;
/// Release constant `NAIF_RHEA`.
pub const NAIF_RHEA: i32 = 605;
/// Release constant `NAIF_TITAN`.
pub const NAIF_TITAN: i32 = 606;
/// Release constant `NAIF_HYPERION`.
pub const NAIF_HYPERION: i32 = 607;
/// Release constant `NAIF_IAPETUS`.
pub const NAIF_IAPETUS: i32 = 608;
/// Release constant `NAIF_MIRANDA`.
pub const NAIF_MIRANDA: i32 = 705;
/// Release constant `NAIF_ARIEL`.
pub const NAIF_ARIEL: i32 = 701;
/// Release constant `NAIF_UMBRIEL`.
pub const NAIF_UMBRIEL: i32 = 702;
/// Release constant `NAIF_TITANIA`.
pub const NAIF_TITANIA: i32 = 703;
/// Release constant `NAIF_OBERON`.
pub const NAIF_OBERON: i32 = 704;
/// Release constant `NAIF_TRITON`.
pub const NAIF_TRITON: i32 = 801;
/// Release constant `NAIF_PHOBOS`.
pub const NAIF_PHOBOS: i32 = 401;
/// Release constant `NAIF_DEIMOS`.
pub const NAIF_DEIMOS: i32 = 402;
/// Release constant `NAIF_CHARON`.
pub const NAIF_CHARON: i32 = 901;
/// Release constant `ASC`.
pub const ASC: i32 = 0;
/// Release constant `MC`.
pub const MC: i32 = 1;
/// Release constant `ARMC`.
pub const ARMC: i32 = 2;
/// Release constant `VERTEX`.
pub const VERTEX: i32 = 3;
/// Release constant `EQUASC`.
pub const EQUASC: i32 = 4;
/// Release constant `COASC1`.
pub const COASC1: i32 = 5;
/// Release constant `COASC2`.
pub const COASC2: i32 = 6;
/// Release constant `POLASC`.
pub const POLASC: i32 = 7;
/// Release constant `NASCMC`.
pub const NASCMC: i32 = 8;
/// Release constant `FIXSTAR`.
pub const FIXSTAR: i32 = -10;
/// Release constant `HARRINGTON`.
pub const HARRINGTON: i32 = 50;
/// Release constant `NEPTUNE_ADAMS`.
pub const NEPTUNE_ADAMS: i32 = 52;
/// Release constant `NEPTUNE_LEVERRIER`.
pub const NEPTUNE_LEVERRIER: i32 = 51;
/// Release constant `NIBIRU`.
pub const NIBIRU: i32 = 49;
/// Release constant `PLUTO_LOWELL`.
pub const PLUTO_LOWELL: i32 = 53;
/// Release constant `PLUTO_PICKERING`.
pub const PLUTO_PICKERING: i32 = 54;
/// Release constant `FICT_MAX`.
pub const FICT_MAX: i32 = 999;
/// Release constant `FICT_OFFSET_1`.
pub const FICT_OFFSET_1: i32 = 39;
/// Release constant `MAX_STNAME`.
pub const MAX_STNAME: i32 = 256;
/// Release constant `NALL_NAT_POINTS`.
pub const NALL_NAT_POINTS: i32 = 38;
/// Release constant `PLMOON_OFFSET`.
pub const PLMOON_OFFSET: i32 = 9000;
/// Release constant `PLMOON_IO`.
pub const PLMOON_IO: i32 = (PLMOON_OFFSET + NAIF_IO);
/// Release constant `PLMOON_EUROPA`.
pub const PLMOON_EUROPA: i32 = (PLMOON_OFFSET + NAIF_EUROPA);
/// Release constant `PLMOON_GANYMEDE`.
pub const PLMOON_GANYMEDE: i32 = (PLMOON_OFFSET + NAIF_GANYMEDE);
/// Release constant `PLMOON_CALLISTO`.
pub const PLMOON_CALLISTO: i32 = (PLMOON_OFFSET + NAIF_CALLISTO);
/// Release constant `PLMOON_MIMAS`.
pub const PLMOON_MIMAS: i32 = (PLMOON_OFFSET + NAIF_MIMAS);
/// Release constant `PLMOON_ENCELADUS`.
pub const PLMOON_ENCELADUS: i32 = (PLMOON_OFFSET + NAIF_ENCELADUS);
/// Release constant `PLMOON_TETHYS`.
pub const PLMOON_TETHYS: i32 = (PLMOON_OFFSET + NAIF_TETHYS);
/// Release constant `PLMOON_DIONE`.
pub const PLMOON_DIONE: i32 = (PLMOON_OFFSET + NAIF_DIONE);
/// Release constant `PLMOON_RHEA`.
pub const PLMOON_RHEA: i32 = (PLMOON_OFFSET + NAIF_RHEA);
/// Release constant `PLMOON_TITAN`.
pub const PLMOON_TITAN: i32 = (PLMOON_OFFSET + NAIF_TITAN);
/// Release constant `PLMOON_HYPERION`.
pub const PLMOON_HYPERION: i32 = (PLMOON_OFFSET + NAIF_HYPERION);
/// Release constant `PLMOON_IAPETUS`.
pub const PLMOON_IAPETUS: i32 = (PLMOON_OFFSET + NAIF_IAPETUS);
/// Release constant `PLMOON_MIRANDA`.
pub const PLMOON_MIRANDA: i32 = (PLMOON_OFFSET + NAIF_MIRANDA);
/// Release constant `PLMOON_ARIEL`.
pub const PLMOON_ARIEL: i32 = (PLMOON_OFFSET + NAIF_ARIEL);
/// Release constant `PLMOON_UMBRIEL`.
pub const PLMOON_UMBRIEL: i32 = (PLMOON_OFFSET + NAIF_UMBRIEL);
/// Release constant `PLMOON_TITANIA`.
pub const PLMOON_TITANIA: i32 = (PLMOON_OFFSET + NAIF_TITANIA);
/// Release constant `PLMOON_OBERON`.
pub const PLMOON_OBERON: i32 = (PLMOON_OFFSET + NAIF_OBERON);
/// Release constant `PLMOON_TRITON`.
pub const PLMOON_TRITON: i32 = (PLMOON_OFFSET + NAIF_TRITON);
/// Release constant `PLMOON_PHOBOS`.
pub const PLMOON_PHOBOS: i32 = (PLMOON_OFFSET + NAIF_PHOBOS);
/// Release constant `PLMOON_DEIMOS`.
pub const PLMOON_DEIMOS: i32 = (PLMOON_OFFSET + NAIF_DEIMOS);
/// Release constant `PLMOON_CHARON`.
pub const PLMOON_CHARON: i32 = (PLMOON_OFFSET + NAIF_CHARON);
/// Release constant `FLG_CENTER_BODY`.
pub const FLG_CENTER_BODY: i32 = 1048576;
/// Release constant `FLG_DEFAULTEPH`.
pub const FLG_DEFAULTEPH: i32 = 2;
/// Release constant `FLG_DPSIDEPS_1980`.
pub const FLG_DPSIDEPS_1980: i32 = 262144;
/// Release constant `FLG_JPLHOR`.
pub const FLG_JPLHOR: i32 = 262144;
/// Release constant `FLG_JPLHOR_APPROX`.
pub const FLG_JPLHOR_APPROX: i32 = 524288;
/// Release constant `FLG_ORBEL_AA`.
pub const FLG_ORBEL_AA: i32 = 32768;
/// Release constant `FLG_TEST_PLMOON`.
pub const FLG_TEST_PLMOON: i32 = 2228280;
/// Release constant `FLG_TROPICAL`.
pub const FLG_TROPICAL: i32 = 0;
/// Release constant `FNAME_DE200`.
pub const FNAME_DE200: &str = "de200.eph";
/// Release constant `FNAME_DE403`.
pub const FNAME_DE403: &str = "de403.eph";
/// Release constant `FNAME_DE404`.
pub const FNAME_DE404: &str = "de404.eph";
/// Release constant `FNAME_DE405`.
pub const FNAME_DE405: &str = "de405.eph";
/// Release constant `FNAME_DE406`.
pub const FNAME_DE406: &str = "de406.eph";
/// Release constant `FNAME_DFT`.
pub const FNAME_DFT: &str = "de431.eph";
/// Release constant `FNAME_DFT2`.
pub const FNAME_DFT2: &str = "de406.eph";
/// Release constant `ECL_HYBRID`.
pub const ECL_HYBRID: i32 = 32;
/// Release constant `ECL_OCC_BEG_DAYLIGHT`.
pub const ECL_OCC_BEG_DAYLIGHT: i32 = 8192;
/// Release constant `ECL_OCC_END_DAYLIGHT`.
pub const ECL_OCC_END_DAYLIGHT: i32 = 16384;
/// Release constant `ECL_PARTBEG_VISIBLE`.
pub const ECL_PARTBEG_VISIBLE: i32 = 512;
/// Release constant `ECL_PARTEND_VISIBLE`.
pub const ECL_PARTEND_VISIBLE: i32 = 4096;
/// Release constant `ECL_PENUMBBEG_VISIBLE`.
pub const ECL_PENUMBBEG_VISIBLE: i32 = 8192;
/// Release constant `ECL_PENUMBEND_VISIBLE`.
pub const ECL_PENUMBEND_VISIBLE: i32 = 16384;
/// Release constant `ECL_TOTBEG_VISIBLE`.
pub const ECL_TOTBEG_VISIBLE: i32 = 1024;
/// Release constant `ECL_TOTEND_VISIBLE`.
pub const ECL_TOTEND_VISIBLE: i32 = 2048;
/// Release constant `NSIDM_PREDEF`.
pub const NSIDM_PREDEF: i32 = 47;
/// Release constant `SIDBITS`.
pub const SIDBITS: i32 = 256;
/// Release constant `SIDBIT_ECL_DATE`.
pub const SIDBIT_ECL_DATE: i32 = 2048;
/// Release constant `SIDBIT_ECL_T0`.
pub const SIDBIT_ECL_T0: i32 = 256;
/// Release constant `SIDBIT_NO_PREC_OFFSET`.
pub const SIDBIT_NO_PREC_OFFSET: i32 = 4096;
/// Release constant `SIDBIT_PREC_ORIG`.
pub const SIDBIT_PREC_ORIG: i32 = 8192;
/// Release constant `SIDBIT_SSY_PLANE`.
pub const SIDBIT_SSY_PLANE: i32 = 512;
/// Release constant `SIDBIT_USER_UT`.
pub const SIDBIT_USER_UT: i32 = 1024;
/// Release constant `MODEL_BIAS`.
pub const MODEL_BIAS: i32 = 4;
/// Release constant `MODEL_DELTAT`.
pub const MODEL_DELTAT: i32 = 0;
/// Release constant `MODEL_JPLHORA_MODE`.
pub const MODEL_JPLHORA_MODE: i32 = 6;
/// Release constant `MODEL_JPLHOR_MODE`.
pub const MODEL_JPLHOR_MODE: i32 = 5;
/// Release constant `MODEL_NUT`.
pub const MODEL_NUT: i32 = 3;
/// Release constant `MODEL_PREC_LONGTERM`.
pub const MODEL_PREC_LONGTERM: i32 = 1;
/// Release constant `MODEL_PREC_SHORTTERM`.
pub const MODEL_PREC_SHORTTERM: i32 = 2;
/// Release constant `MODEL_SIDT`.
pub const MODEL_SIDT: i32 = 7;
/// Release constant `NSE_MODELS`.
pub const NSE_MODELS: i32 = 8;
/// Release constant `MOD_BIAS_DEFAULT`.
pub const MOD_BIAS_DEFAULT: i32 = 3;
/// Release constant `MOD_BIAS_IAU2000`.
pub const MOD_BIAS_IAU2000: i32 = 2;
/// Release constant `MOD_BIAS_IAU2006`.
pub const MOD_BIAS_IAU2006: i32 = 3;
/// Release constant `MOD_BIAS_NONE`.
pub const MOD_BIAS_NONE: i32 = 1;
/// Release constant `MOD_NBIAS`.
pub const MOD_NBIAS: i32 = 3;
/// Release constant `MOD_DELTAT_DEFAULT`.
pub const MOD_DELTAT_DEFAULT: i32 = 5;
/// Release constant `MOD_DELTAT_ESPENAK_MEEUS_2006`.
pub const MOD_DELTAT_ESPENAK_MEEUS_2006: i32 = 4;
/// Release constant `MOD_DELTAT_STEPHENSON_1997`.
pub const MOD_DELTAT_STEPHENSON_1997: i32 = 2;
/// Release constant `MOD_DELTAT_STEPHENSON_ETC_2016`.
pub const MOD_DELTAT_STEPHENSON_ETC_2016: i32 = 5;
/// Release constant `MOD_DELTAT_STEPHENSON_MORRISON_1984`.
pub const MOD_DELTAT_STEPHENSON_MORRISON_1984: i32 = 1;
/// Release constant `MOD_DELTAT_STEPHENSON_MORRISON_2004`.
pub const MOD_DELTAT_STEPHENSON_MORRISON_2004: i32 = 3;
/// Release constant `MOD_NDELTAT`.
pub const MOD_NDELTAT: i32 = 5;
/// Release constant `MOD_JPLHORA_1`.
pub const MOD_JPLHORA_1: i32 = 1;
/// Release constant `MOD_JPLHORA_2`.
pub const MOD_JPLHORA_2: i32 = 2;
/// Release constant `MOD_JPLHORA_3`.
pub const MOD_JPLHORA_3: i32 = 3;
/// Release constant `MOD_JPLHORA_DEFAULT`.
pub const MOD_JPLHORA_DEFAULT: i32 = 3;
/// Release constant `MOD_NJPLHORA`.
pub const MOD_NJPLHORA: i32 = 3;
/// Release constant `MOD_JPLHOR_DEFAULT`.
pub const MOD_JPLHOR_DEFAULT: i32 = 1;
/// Release constant `MOD_JPLHOR_LONG_AGREEMENT`.
pub const MOD_JPLHOR_LONG_AGREEMENT: i32 = 1;
/// Release constant `MOD_NJPLHOR`.
pub const MOD_NJPLHOR: i32 = 2;
/// Release constant `MOD_NUT_DEFAULT`.
pub const MOD_NUT_DEFAULT: i32 = 4;
/// Release constant `MOD_NUT_IAU_1980`.
pub const MOD_NUT_IAU_1980: i32 = 1;
/// Release constant `MOD_NUT_IAU_2000A`.
pub const MOD_NUT_IAU_2000A: i32 = 3;
/// Release constant `MOD_NUT_IAU_2000B`.
pub const MOD_NUT_IAU_2000B: i32 = 4;
/// Release constant `MOD_NUT_IAU_CORR_1987`.
pub const MOD_NUT_IAU_CORR_1987: i32 = 2;
/// Release constant `MOD_NUT_WOOLARD`.
pub const MOD_NUT_WOOLARD: i32 = 5;
/// Release constant `MOD_NNUT`.
pub const MOD_NNUT: i32 = 5;
/// Release constant `MOD_PREC_BRETAGNON_2003`.
pub const MOD_PREC_BRETAGNON_2003: i32 = 7;
/// Release constant `MOD_PREC_DEFAULT`.
pub const MOD_PREC_DEFAULT: i32 = 9;
/// Release constant `MOD_PREC_DEFAULT_SHORT`.
pub const MOD_PREC_DEFAULT_SHORT: i32 = 9;
/// Release constant `MOD_PREC_IAU_1976`.
pub const MOD_PREC_IAU_1976: i32 = 1;
/// Release constant `MOD_PREC_IAU_2000`.
pub const MOD_PREC_IAU_2000: i32 = 6;
/// Release constant `MOD_PREC_IAU_2006`.
pub const MOD_PREC_IAU_2006: i32 = 8;
/// Release constant `MOD_PREC_LASKAR_1986`.
pub const MOD_PREC_LASKAR_1986: i32 = 2;
/// Release constant `MOD_PREC_NEWCOMB`.
pub const MOD_PREC_NEWCOMB: i32 = 11;
/// Release constant `MOD_PREC_OWEN_1990`.
pub const MOD_PREC_OWEN_1990: i32 = 10;
/// Release constant `MOD_PREC_SIMON_1994`.
pub const MOD_PREC_SIMON_1994: i32 = 5;
/// Release constant `MOD_PREC_VONDRAK_2011`.
pub const MOD_PREC_VONDRAK_2011: i32 = 9;
/// Release constant `MOD_PREC_WILLIAMS_1994`.
pub const MOD_PREC_WILLIAMS_1994: i32 = 4;
/// Release constant `MOD_PREC_WILL_EPS_LASK`.
pub const MOD_PREC_WILL_EPS_LASK: i32 = 3;
/// Release constant `MOD_NPREC`.
pub const MOD_NPREC: i32 = 11;
/// Release constant `BIT_FORCE_SLOW_METHOD`.
pub const BIT_FORCE_SLOW_METHOD: i32 = 32768;
/// Release constant `BIT_GEOCTR_NO_ECL_LAT`.
pub const BIT_GEOCTR_NO_ECL_LAT: i32 = 128;
/// Release constant `BIT_HINDU_RISING`.
pub const BIT_HINDU_RISING: i32 = 896;
/// Release constant `ACRONYCHAL_RISING`.
pub const ACRONYCHAL_RISING: i32 = 5;
/// Release constant `ACRONYCHAL_SETTING`.
pub const ACRONYCHAL_SETTING: i32 = 6;
/// Release constant `COSMICAL_SETTING`.
pub const COSMICAL_SETTING: i32 = 6;
/// Release constant `TIDAL_26` (bit-exact).
pub const TIDAL_26: f64 = f64::from_bits(0xC03A000000000000);
/// Release constant `TIDAL_JPLEPH` (bit-exact).
pub const TIDAL_JPLEPH: f64 = f64::from_bits(0xC039CCCCCCCCCCCD);
/// Release constant `TIDAL_MOSEPH` (bit-exact).
pub const TIDAL_MOSEPH: f64 = f64::from_bits(0xC039947AE147AE14);
/// Release constant `TIDAL_STEPHENSON_2016` (bit-exact).
pub const TIDAL_STEPHENSON_2016: f64 = f64::from_bits(0xC039D9999999999A);
/// Release constant `TIDAL_SWIEPH` (bit-exact).
pub const TIDAL_SWIEPH: f64 = f64::from_bits(0xC039CCCCCCCCCCCD);
/// Release constant `ASTNAMFILE`.
pub const ASTNAMFILE: &str = "seasnam.txt";
/// Release constant `EPHE_PATH`.
pub const EPHE_PATH: &str = ".:/users/ephe2/:/users/ephe/";
/// Release constant `FICTFILE`.
pub const FICTFILE: &str = "seorbel.txt";
/// Release constant `FNAME_DE431`.
pub const FNAME_DE431: &str = "de431.eph";
/// Release constant `SE_FNAME_DE431`: intentional alias of `FNAME_DE431`,
/// matching the release assignment (asserted in tests).
pub const SE_FNAME_DE431: &str = "de431.eph";
/// Release constant `STARFILE`.
pub const STARFILE: &str = "sefstars.txt";
/// Release constant `STARFILE_OLD`.
pub const STARFILE_OLD: &str = "fixstars.cat";
/// Release constant `AUNIT_TO_KM` (bit-exact).
pub const AUNIT_TO_KM: f64 = f64::from_bits(0x41A1D55D5D666666);
/// Release constant `AUNIT_TO_LIGHTYEAR` (bit-exact).
pub const AUNIT_TO_LIGHTYEAR: f64 = f64::from_bits(0x3EF094A33C2EA847);
/// Release constant `AUNIT_TO_PARSEC` (bit-exact).
pub const AUNIT_TO_PARSEC: f64 = f64::from_bits(0x3ED455A5B2FF8F38);
/// Release constant `DE_NUMBER`.
pub const DE_NUMBER: i32 = 431;
/// Release constant `DELTAT_AUTOMATIC` (bit-exact).
pub const DELTAT_AUTOMATIC: f64 = f64::from_bits(0xBDDB7CDFD9D7BDBB);
/// Release constant `MIXEDOPIC_FLAG`.
pub const MIXEDOPIC_FLAG: i32 = 2;
/// Release constant `PHOTOPIC_FLAG`.
pub const PHOTOPIC_FLAG: i32 = 0;
/// Release constant `SCOTOPIC_FLAG`.
pub const SCOTOPIC_FLAG: i32 = 1;
/// Release constant `SIMULATE_VICTORVB`.
pub const SIMULATE_VICTORVB: i32 = 1;
/// Release constant `TJD_INVALID` (bit-exact).
pub const TJD_INVALID: f64 = f64::from_bits(0x4197D783FC000000);
/// Horizons/NAIF identity map (key, Horizons id, NAIF id), sorted by key.
/// Shape mirrors the release dict; unmapped ids yield `None` via the lookup fns.
pub const SPK_BODY_NAME_MAP: &[(i32, &str, i32)] = &[
    (15, "2060", 2002060),
    (16, "5145", 2005145),
    (17, "Ceres;", 2000001),
    (18, "Pallas;", 2000002),
    (19, "Juno;", 2000003),
    (20, "Vesta;", 2000004),
    (10010, "Hygiea;", 2000010),
    (10016, "16", 2000016),
    (10052, "52", 2000052),
    (10055, "55", 2000055),
    (10080, "80", 2000080),
    (10087, "87", 2000087),
    (10433, "433", 2000433),
    (10511, "Davida;", 2000511),
    (10704, "Interamnia;", 2000704),
    (10944, "944", 2000944),
    (11181, "1181", 2001181),
    (11221, "1221", 2001221),
    (11566, "1566", 2001566),
    (11685, "1685", 2001685),
    (14179, "4179", 2004179),
    (17066, "7066", 2007066),
    (18405, "8405", 2008405),
    (20199, "10199", 2010199),
    (30000, "20000", 2020000),
    (35143, "25143", 2025143),
    (38978, "28978", 2028978),
    (60000, "50000", 2050000),
    (100377, "90377", 2090377),
    (100482, "90482", 2090482),
    (109942, "99942", 2099942),
    (111955, "Bennu;", 2101955),
    (146108, "136108", 2136108),
    (146199, "136199", 2136199),
    (146472, "136472", 2136472),
    (172173, "162173", 2162173),
    (235088, "225088", 2225088),
];
/// Bodies blocked from automatic SPK download, with reason.
/// Blocked bodies stay mapped for manual kernel registration; auto-download path must skip them.
pub const SPK_AUTO_DOWNLOAD_BLOCKED: &[(i32, &str)] = &[(
    111955,
    "JPL Horizons blocks SPK generation for Bennu. Use ASSIST n-body or register a mission-specific kernel (e.g. bennu_refdrmc_v1.bsp from OSIRIS-REx).",
)];
/// Bodies required for high-precision calculations, sorted.
/// Every entry is mapped in `SPK_BODY_NAME_MAP`; used to drive SPK ensure/download.
pub const REQUIRED_SPK_BODIES: &[i32] = &[15, 16, 17, 18, 19, 20, 17066, 146199];

/// JPL Horizons target identifier for an integer libephemeris body id (`ipl`); `None` when unmapped.
///
/// Mirrors release `get_horizons_id`: the Horizons string from
/// [`SPK_BODY_NAME_MAP`], or `None` when the body is not in the mapping.
pub fn get_horizons_id(ipl: i32) -> Option<&'static str> {
    get_spk_body_info_from_map(ipl).map(|(horizons_id, _)| horizons_id)
}

/// NAIF SPICE id for an integer libephemeris body id (`ipl`); `None` when unmapped.
///
/// Mirrors release `get_naif_id_from_ipl`: the NAIF id from
/// [`SPK_BODY_NAME_MAP`], or `None` when the body is not in the mapping.
pub fn get_naif_id_from_ipl(ipl: i32) -> Option<i32> {
    get_spk_body_info_from_map(ipl).map(|(_, naif_id)| naif_id)
}

/// Horizons id and NAIF id for an integer libephemeris body id (`ipl`); `None` when unmapped.
///
/// Mirrors release `get_spk_body_info_from_map`: the `(horizons_id, naif_id)`
/// pair from [`SPK_BODY_NAME_MAP`], or `None` when the body is not in the
/// mapping. [`SPK_BODY_NAME_MAP`] is sorted by key, so the lookup is an exact
/// binary search with the same membership semantics as the release dict.
pub fn get_spk_body_info_from_map(ipl: i32) -> Option<(&'static str, i32)> {
    SPK_BODY_NAME_MAP
        .binary_search_by_key(&ipl, |(key, _, _)| *key)
        .ok()
        .map(|index| {
            let (_, horizons_id, naif_id) = SPK_BODY_NAME_MAP[index];
            (horizons_id, naif_id)
        })
}

/// Whether JPL Horizons blocks SPK generation for an integer body id (`ipl`); `false` when unmapped.
///
/// Mirrors release `is_spk_auto_download_blocked`: membership in
/// [`SPK_AUTO_DOWNLOAD_BLOCKED`]. Bodies in that table keep valid
/// Horizons/NAIF identifiers for manual kernel registration but must skip
/// the automatic SPK download path.
pub fn is_spk_auto_download_blocked(ipl: i32) -> bool {
    SPK_AUTO_DOWNLOAD_BLOCKED
        .binary_search_by_key(&ipl, |(key, _)| *key)
        .is_ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn release_table_count_matches() {
        // Every public release name is implemented exactly once.
        assert_eq!(RELEASE_CONST_COUNT, 594);
    }

    #[test]
    fn spk_lookup_functions_match_release_contract() {
        // Pinned release behaviors (verified against constants.py @ 37ffb05):
        // get_horizons_id(CHIRON) == '2060', get_naif_id_from_ipl(ERIS) ==
        // 2136199, get_spk_body_info_from_map(ERIS) == ('136199', 2136199),
        // SUN and unknown ids map to None, BENNU is download-blocked.
        assert_eq!(get_horizons_id(CHIRON), Some("2060"));
        assert_eq!(get_horizons_id(SUN), None);
        assert_eq!(get_horizons_id(999_999), None);
        assert_eq!(get_naif_id_from_ipl(CHIRON), Some(2002060));
        assert_eq!(get_naif_id_from_ipl(ERIS), Some(2136199));
        assert_eq!(get_naif_id_from_ipl(SUN), None);
        assert_eq!(get_spk_body_info_from_map(ERIS), Some(("136199", 2136199)));
        assert_eq!(get_spk_body_info_from_map(SUN), None);
        assert!(is_spk_auto_download_blocked(BENNU));
        assert!(!is_spk_auto_download_blocked(CHIRON));
        assert!(!is_spk_auto_download_blocked(999_999));
    }

    #[test]
    fn float_constants_are_bit_exact() {
        // Representative bit identities, verified with pure-Python
        // struct packing of the release decimal literals (no package
        // import): every Rust `from_bits` word below matched first run.
        assert_eq!(AUNIT.to_bits(), 0x41A1D55D5D666666);
        assert_eq!(AUNIT_TO_KM.to_bits(), 0x41A1D55D5D666666);
        assert_eq!(TIDAL_DE200.to_bits(), 0xC037E504816F0069);
        assert_eq!(TIDAL_DE403.to_bits(), 0xC039947AE147AE14);
        assert_eq!(TIDAL_DE405.to_bits(), 0xC039D374BC6A7EFA);
        assert_eq!(TIDAL_DE421.to_bits(), 0xC039D9999999999A);
        assert_eq!(TIDAL_DE430.to_bits(), 0xC039D1EB851EB852);
        assert_eq!(TIDAL_DE431.to_bits(), 0xC039CCCCCCCCCCCD);
        assert_eq!(TIDAL_DEFAULT.to_bits(), 0xC039CCCCCCCCCCCD);
        assert_eq!(TIDAL_DE440.to_bits(), 0xC039EF9DB22D0E56);
        assert_eq!(J2000.to_bits(), 0x4142B42C80000000);
        assert_eq!(J1991_25.to_bits(), 0x4142ADEE88000000);
        assert_eq!(J1900.to_bits(), 0x41426CD600000000);
        assert_eq!(B1950.to_bits(), 0x414290813633E7F9);
        assert_eq!(DAYS_PER_JULIAN_YEAR.to_bits(), 0x4076D40000000000);
        assert_eq!(DAYS_PER_JULIAN_CENTURY.to_bits(), 0x40E1D5A000000000);
        assert_eq!(TIDAL_26.to_bits(), 0xC03A000000000000);
        assert_eq!(TIDAL_JPLEPH.to_bits(), 0xC039CCCCCCCCCCCD);
        assert_eq!(TIDAL_MOSEPH.to_bits(), 0xC039947AE147AE14);
        assert_eq!(TIDAL_STEPHENSON_2016.to_bits(), 0xC039D9999999999A);
        assert_eq!(AUNIT_TO_LIGHTYEAR.to_bits(), 0x3EF094A33C2EA847);
        assert_eq!(AUNIT_TO_PARSEC.to_bits(), 0x3ED455A5B2FF8F38);
        assert_eq!(DELTAT_AUTOMATIC.to_bits(), 0xBDDB7CDFD9D7BDBB);
        assert_eq!(TJD_INVALID.to_bits(), 0x4197D783FC000000);
        // Aliased floats share the exact word by construction.
        assert_eq!(TIDAL_DE404.to_bits(), TIDAL_DE403.to_bits());
        assert_eq!(TIDAL_DE406.to_bits(), TIDAL_DE405.to_bits());
        assert_eq!(TIDAL_DE422.to_bits(), TIDAL_DE421.to_bits());
        assert_eq!(TIDAL_SWIEPH.to_bits(), TIDAL_JPLEPH.to_bits());
    }

    #[test]
    fn integer_expressions_strings_and_compat_names_match() {
        // Expression-preserved offsets evaluate to the release values.
        assert_eq!(VARUNA, 30_000);
        assert_eq!(CUPIDO, 40);
        assert_eq!(ERIS, 146_199);
        assert_eq!(REGULUS, 1_000_001);
        assert_eq!(FLG_ASTROMETRIC, FLG_NOABERR | FLG_NOGDEFL);
        assert_eq!(
            ECL_ALLTYPES_SOLAR,
            ECL_CENTRAL
                | ECL_NONCENTRAL
                | ECL_TOTAL
                | ECL_ANNULAR
                | ECL_PARTIAL
                | ECL_ANNULAR_TOTAL
        );
        assert_eq!(BIT_HINDU_RISING, 896);
        // Aliases without a separate numeric identity.
        assert_eq!(TRANSPLUTO, ISIS);
        assert_eq!(SELENA, WHITE_MOON);
        assert_eq!(MORNING_FIRST, HELIACAL_RISING);
        assert_eq!(COSMICAL_SETTING, ACRONYCHAL_SETTING);
        assert_eq!(TIDAL_DEFAULT, TIDAL_DE431);
        assert_eq!(PLMOON_OFFSET, MOON_OFFSET);
        assert_eq!(PLMOON_IO, PLMOON_OFFSET + NAIF_IO);
        // String/file policy is, including the single SE_ exception.
        assert_eq!(FNAME_DE431, "de431.eph");
        assert_eq!(SE_FNAME_DE431, FNAME_DE431);
        assert_eq!(EPHE_PATH, ".:/users/ephe2/:/users/ephe/");
        assert_eq!(FICTFILE, "seorbel.txt");
        // Sorted-table invariants for binary-search lookups.
        assert!(SPK_BODY_NAME_MAP.windows(2).all(|w| w[0].0 < w[1].0));
        assert!(REQUIRED_SPK_BODIES.windows(2).all(|w| w[0] < w[1]));
    }

    #[test]
    fn spk_lookup_functions_agree_on_every_mapped_body() {
        // The three lookup functions share one table, so they must agree on
        // every mapped key, and every required/blocked body must be mapped
        // (or explicitly unmapped, which the release forbids by construction
        // of REQUIRED_SPK_BODIES from mapped bodies).
        for &(key, horizons_id, naif_id) in SPK_BODY_NAME_MAP {
            assert_eq!(get_horizons_id(key), Some(horizons_id));
            assert_eq!(get_naif_id_from_ipl(key), Some(naif_id));
            assert_eq!(
                get_spk_body_info_from_map(key),
                Some((horizons_id, naif_id))
            );
        }
        for &body in REQUIRED_SPK_BODIES {
            assert!(
                get_spk_body_info_from_map(body).is_some(),
                "required body {body} must be mapped"
            );
        }
        for &(body, _) in SPK_AUTO_DOWNLOAD_BLOCKED {
            assert!(
                get_spk_body_info_from_map(body).is_some(),
                "blocked body {body} must stay mapped for manual registration"
            );
        }
    }
}
