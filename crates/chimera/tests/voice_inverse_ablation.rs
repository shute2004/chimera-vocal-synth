//! Diagnostic-only ablation of glottal area, coupled port flow, and radiated PCM.
//! CHIMERA synthesis is unchanged; harmonic ratios are not additive under coupling.
use chimera::airway::frequency::transfer_function;
use chimera::anatomy::human::{ArticKey, BodyMix};
use chimera::phonetics::planner::vowel_artic;
use chimera::voice::human_voice::{HumanVoice, SourcePatch};

const FS: f64 = 48_000.0;
const PITCH: f64 = 420.0;
const N: usize = 8192;

#[derive(Clone, Copy)]
struct Settings {
    female: f64,
    jaw: f64,
    tongue_y: f64,
    tongue_x: f64,
    bright: f64,
    tension: f64,
    pressure: f64,
    breath: f64,
}
impl Settings {
    fn baseline() -> Self {
        // Generic reproducible setting; not a fitted human recording.
        Self {
            female: 0.35,
            jaw: 0.15,
            tongue_y: 0.0,
            tongue_x: 0.0,
            bright: 0.75,
            tension: 0.40,
            pressure: 8.0,
            breath: 0.08,
        }
    }
}
struct Measure {
    area: f64,
    flow: f64,
    output: f64,
    tract: f64,
    output_h1h2: f64,
}
fn harmonic(samples: &[f64], n: usize) -> f64 {
    let samples = &samples[samples.len() - N..];
    let mut re = 0.0;
    let mut im = 0.0;
    for (i, &v) in samples.iter().enumerate() {
        let t = i as f64;
        let window = 0.5 - 0.5 * (std::f64::consts::TAU * t / (N - 1) as f64).cos();
        let phase = std::f64::consts::TAU * (PITCH * n as f64 / FS) * t;
        re += v * window * phase.cos();
        im -= v * window * phase.sin();
    }
    re * re + im * im
}
fn ratio(samples: &[f64], a: usize, b: usize) -> f64 {
    10.0 * ((harmonic(samples, a) + 1e-20) / (harmonic(samples, b) + 1e-20)).log10()
}
fn measure(s: Settings) -> Measure {
    let body = BodyMix::new(s.female, 0.0);
    let mut voice = HumanVoice::new(FS, body.body(), 2);
    let mut artic = vowel_artic("a", body).unwrap();
    artic.set(ArticKey::Jaw, artic.get(ArticKey::Jaw).unwrap() + s.jaw);
    artic.set(
        ArticKey::TbY,
        artic.get(ArticKey::TbY).unwrap() + s.tongue_y,
    );
    artic.set(
        ArticKey::TbX,
        artic.get(ArticKey::TbX).unwrap() + s.tongue_x,
    );
    voice.set_articulation(&artic, None);
    voice.set_source(&SourcePatch {
        f0: Some(PITCH),
        glide: Some(0.025),
        brightness: Some(s.bright),
        tension: Some(s.tension),
        pressure: Some(s.pressure),
        breathiness: Some(s.breath),
        roughness: Some(0.0),
        drift: Some(0.0),
        vibrato_depth: Some(0.0),
        ..Default::default()
    });
    voice.set_gate(1.0);
    let mut output = Vec::with_capacity(20_000);
    let mut glottal_area = Vec::with_capacity(20_000);
    let mut glottal_flow = Vec::with_capacity(20_000);
    for _ in 0..20_000 {
        let mut sample = [0.0f32; 1];
        voice.render(&mut sample);
        output.push(sample[0] as f64);
        glottal_area.push(voice.glottis.area);
        glottal_flow.push(voice.airway.net.u[voice.airway.glottis_port]);
    }
    let tf = transfer_function(
        &voice.airway.net.linear_elements(),
        voice.airway.pharynx.first_node(),
        &[PITCH * 2.0, PITCH * 4.0],
        &[voice.airway.glottis_port],
    );
    let tract = 20.0 * (tf.mag[1] / tf.mag[0]).log10();
    Measure {
        area: ratio(&glottal_area, 4, 2),
        flow: ratio(&glottal_flow, 4, 2),
        output: ratio(&output, 4, 2),
        output_h1h2: ratio(&output, 1, 2),
        tract,
    }
}
#[test]
fn source_tract_coupling_ablation() {
    let base = Settings::baseline();
    let b = measure(base);
    let cases = [
        ("jaw_more", Settings { jaw: 0.55, ..base }),
        ("jaw_less", Settings { jaw: -0.40, ..base }),
        (
            "tongue_up",
            Settings {
                tongue_y: 0.20,
                ..base
            },
        ),
        (
            "tongue_down",
            Settings {
                tongue_y: -0.20,
                ..base
            },
        ),
        (
            "soft_closure",
            Settings {
                bright: 0.10,
                ..base
            },
        ),
        (
            "pressed_closure",
            Settings {
                tension: 0.80,
                ..base
            },
        ),
        (
            "pressure_low",
            Settings {
                pressure: 5.0,
                ..base
            },
        ),
        (
            "pressure_high",
            Settings {
                pressure: 11.0,
                ..base
            },
        ),
        (
            "breathy",
            Settings {
                breath: 0.45,
                ..base
            },
        ),
        (
            "jaw_and_closure",
            Settings {
                jaw: 0.55,
                bright: 0.10,
                ..base
            },
        ),
    ];
    println!("case,area_H4_H2,flow_H4_H2,output_H4_H2,linear_tract_H4_H2,output_H1_H2");
    let print = |name: &str, m: &Measure| {
        println!(
            "{name},{:.4},{:.4},{:.4},{:.4},{:.4}",
            m.area, m.flow, m.output, m.tract, m.output_h1h2
        );
        assert!(
            [m.area, m.flow, m.output, m.tract, m.output_h1h2]
                .iter()
                .all(|x| x.is_finite())
        );
    };
    print("baseline", &b);
    let mut results = Vec::new();
    for (name, s) in cases {
        let m = measure(s);
        print(name, &m);
        results.push((name, m));
    }
    let jaw = results
        .iter()
        .find(|(name, _)| *name == "jaw_more")
        .unwrap()
        .1
        .output;
    let soft = results
        .iter()
        .find(|(name, _)| *name == "soft_closure")
        .unwrap()
        .1
        .output;
    let joint = results
        .iter()
        .find(|(name, _)| *name == "jaw_and_closure")
        .unwrap()
        .1
        .output;
    let interaction = joint - jaw - soft + b.output;
    println!("source_x_tract_interaction_H4H2_dB={interaction:.6}");
    assert!(
        interaction.abs() > 0.01,
        "Coupled response must not be assumed additive"
    );
}

/// Run with --ignored --nocapture; scans the physical engine, not a toy filter.
/// The grid is bounded, so failure to hit a desired spectrum is not proof of impossibility.
#[test]
#[ignore = "diagnostic grid: 75 rendered vocal-tract configurations"]
fn scan_tongue_jaw_resonance_range() {
    let base = Settings::baseline();
    let mut rows = Vec::new();
    for jaw in [-0.7, -0.35, 0.0, 0.35, 0.7] {
        for tongue_y in [-0.3, -0.15, 0.0, 0.15, 0.3] {
            for tongue_x in [-0.3, 0.0, 0.3] {
                let measurement = measure(Settings {
                    jaw,
                    tongue_y,
                    tongue_x,
                    ..base
                });
                rows.push((
                    measurement.output,
                    jaw,
                    tongue_y,
                    tongue_x,
                    measurement.flow,
                    measurement.area,
                    measurement.tract,
                ));
            }
        }
    }
    assert_eq!(rows.len(), 75);
    assert!(rows.iter().all(|row| row.0.is_finite()));
    rows.sort_by(|a, b| b.0.total_cmp(&a.0));
    println!("best output H4/H2 settings (dB; jaw/tongueY/tongueX):");
    for row in rows.iter().take(8) {
        println!(
            "output={:.3}, jaw={:.2}, tongueY={:.2}, tongueX={:.2}, glottalFlow={:.3}, glottalArea={:.3}, tractLinear={:.3}",
            row.0, row.1, row.2, row.3, row.4, row.5, row.6
        );
    }
    println!(
        "worst={:.3}dB, best={:.3}dB, configs_above_0dB={}",
        rows.last().unwrap().0,
        rows.first().unwrap().0,
        rows.iter().filter(|r| r.0 > 0.0).count()
    );
    assert!(
        rows[0].0 - rows.last().unwrap().0 > 1.0,
        "Vocal-tract sweep did not move spectrum"
    );
}
