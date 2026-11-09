//! How the human voice's source settings shape the glottis.
use chimera::voice::human_voice::{HumanVoice, SourcePatch};
use chimera::anatomy::human::BodyMix;
use chimera::phonetics::planner::vowel_artic;
use chimera::analysis::pitch::{yin_pitch, YinOpts};

#[test]
fn a_pressed_voice_still_closes_abruptly() {
    // a sung or spoken voice closes softly (tension 0.5-0.55), but the "Pressed" control's top end
    // keeps the abrupt closure of a pressed voice, and the closure sharpens with tension throughout
    let mut v = HumanVoice::with_defaults(48000.0);
    let mut sharpness = |tension: f64| {
        v.set_source(&SourcePatch { tension: Some(tension), ..Default::default() });
        v.glottis.params.closure_sharpness
    };
    assert!(sharpness(0.55) < 0.25);
    assert!(sharpness(1.0) > 0.9);
    let mut last = -1.0;
    for i in 0..=20 {
        let t = i as f64 / 20.0;
        let s = sharpness(t);
        assert!(s > last, "closure sharpness falls at tension {t}");
        last = s;
    }
}

#[test]
fn brightness_tilts_the_voice_and_its_default_changes_nothing() {
    // the Bright control: at its default the closure is what tension alone gives; turned up, the
    // folds close more abruptly and the voice carries more high-frequency energy
    let mut v = HumanVoice::with_defaults(48000.0);
    v.set_source(&SourcePatch { tension: Some(0.55), ..Default::default() });
    let default = v.glottis.params.closure_sharpness;
    assert!((default - 0.22).abs() < 1e-12, "{default}");
    let mut tilt = |b: f64| {
        let mut v = HumanVoice::with_defaults(48000.0);
        v.set_source(&SourcePatch { tension: Some(0.55), brightness: Some(b), ..Default::default() });
        v.set_gate(1.0);
        let mut out = vec![0.0f32; 48000];
        for c in out.chunks_mut(128) {
            v.render(c);
        }
        // mean squared first difference over mean square: rises with the share of high frequencies
        let x = &out[24000..];
        let d: f64 = x.windows(2).map(|w| ((w[1] - w[0]) as f64).powi(2)).sum();
        let e: f64 = x.iter().map(|s| (*s as f64).powi(2)).sum();
        (v.glottis.params.closure_sharpness, 10.0 * (d / e).log10())
    };
    let (dark, mid, bright) = (tilt(0.0), tilt(0.25), tilt(1.0));
    assert_eq!(mid.0, default);
    assert!(dark.0 < mid.0 && mid.0 < bright.0);
    println!("high-frequency share (dB): dark {:.1}, default {:.1}, bright {:.1}", dark.1, mid.1, bright.1);
    // little room to darken (the default closure is already soft), a lot to brighten
    assert!(dark.1 < mid.1 && mid.1 < bright.1 - 3.0, "{dark:?} {mid:?} {bright:?}");
}

#[test]
fn zero_roughness_is_a_clean_periodic_reference() {
    let mut v = HumanVoice::with_defaults(48_000.0);
    for (roughness, expected_jitter, expected_shimmer) in [
        (0.0, 0.0, 0.0),
        (0.05, 0.002275, 0.021),
        (0.1, 0.00455, 0.042),
        (0.5, 0.01435, 0.13),
    ] {
        v.set_source(&SourcePatch {
            roughness: Some(roughness),
            ..Default::default()
        });
        assert!((v.glottis.params.jitter - expected_jitter).abs() < 1e-12);
        assert!((v.glottis.params.shimmer - expected_shimmer).abs() < 1e-12);
    }

    // Verify this is audible in the generated waveform, not just the parameter state.
    let body = BodyMix::new(0.25, 0.0);
    let mut voice = HumanVoice::new(48_000.0, body.body(), 2);
    let artic = vowel_artic("ɑ", body).unwrap();
    voice.set_articulation(&artic, None);
    voice.set_source(&SourcePatch {
        f0: Some(150.0),
        brightness: Some(0.4),
        breathiness: Some(0.08),
        tension: Some(0.55),
        pressure: Some(8.0),
        roughness: Some(0.0),
        drift: Some(0.0),
        vibrato_depth: Some(0.0),
        ..Default::default()
    });
    voice.set_gate(1.0);
    let mut pcm = vec![0.0f32; 48_000];
    for block in pcm.chunks_mut(128) {
        voice.render(block);
    }
    let mut aperiodicity = 0.0;
    for start in [20_000, 25_000, 30_000, 35_000] {
        let frame: Vec<f64> = pcm[start..start + 4096].iter().map(|&x| x as f64).collect();
        let p = yin_pitch(
            &frame,
            48_000.0,
            YinOpts {
                min_f0: Some(75.0),
                max_f0: Some(350.0),
                threshold: None,
            },
        );
        assert!((p.f0 - 150.0).abs() < 1.0);
        aperiodicity += p.aperiodicity / 4.0;
    }
    assert!(
        aperiodicity < 0.001,
        "clean periodic voice has noise: {aperiodicity}"
    );
}
