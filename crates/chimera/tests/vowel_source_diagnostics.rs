//! Reproducible isolated vowel comparison, with no reference recording.
use chimera::airway::frequency::{find_peaks, transfer_function};
use chimera::analysis::fft::{Fft, hann};
use chimera::anatomy::human::BodyMix;
use chimera::phonetics::planner::vowel_artic;
use chimera::voice::human_voice::{HumanVoice, SourcePatch};

const FS: f64 = 48_000.0;
const N: usize = 8192;

fn pcm(ipa: &str, brightness: f64, tension: f64) -> (Vec<f32>, Vec<f64>) {
    let body = BodyMix::new(0.25, 0.0);
    let mut voice = HumanVoice::new(FS, body.body(), 2);
    voice.set_articulation(&vowel_artic(ipa, body).unwrap(), None);
    voice.set_source(&SourcePatch {
        f0: Some(150.0),
        brightness: Some(brightness),
        tension: Some(tension),
        pressure: Some(8.0),
        breathiness: Some(0.08),
        roughness: Some(0.0),
        drift: Some(0.0),
        vibrato_depth: Some(0.0),
        ..Default::default()
    });
    voice.set_gate(1.0);
    let mut samples = vec![0.0f32; 48_000];
    for chunk in samples.chunks_mut(128) {
        voice.render(chunk);
    }
    let frequencies: Vec<f64> = (10..=400).map(|i| i as f64 * 10.0).collect();
    let tf = transfer_function(
        &voice.airway.net.linear_elements(),
        voice.airway.pharynx.first_node(),
        &frequencies,
        &[voice.airway.glottis_port],
    );
    let formants = find_peaks(&frequencies, &tf.mag, 150.0)
        .into_iter()
        .take(3)
        .map(|p| p.freq)
        .collect();
    (samples, formants)
}

fn spectrum(sample: &[f32]) -> ([f64; 5], f64) {
    let fft = Fft::new(N);
    let window = hann(N);
    let mut real: Vec<f64> = sample[30_000..30_000 + N]
        .iter()
        .zip(window)
        .map(|(&x, w)| x as f64 * w)
        .collect();
    let mut imag = vec![0.0; N];
    fft.transform(&mut real, &mut imag);
    let power: Vec<f64> = (0..N / 2)
        .map(|k| real[k] * real[k] + imag[k] * imag[k])
        .collect();
    let total = power[1..((8000.0 * N as f64 / FS) as usize)]
        .iter()
        .sum::<f64>();
    let ranges = [
        (0.0, 500.0),
        (500.0, 1000.0),
        (1000.0, 2000.0),
        (2000.0, 4000.0),
        (4000.0, 8000.0),
    ];
    let band = ranges.map(|(a, b)| {
        let lo = ((a * N as f64 / FS).ceil() as usize).max(1);
        let hi = (b * N as f64 / FS).ceil() as usize;
        10.0 * (power[lo..hi].iter().sum::<f64>() / total).log10()
    });
    let harmonics = [150.0, 300.0, 450.0].map(|hz| {
        let k = (hz * N as f64 / FS).round() as usize;
        power[k - 1..=k + 1].iter().sum::<f64>()
    });
    let h1h2 = 10.0 * (harmonics[0] / harmonics[1]).log10();

    (band, h1h2)
}

#[test]
fn initial_vowel_and_source_measurement() {
    let mut outcomes = Vec::new();
    for (vowel, bright, tension) in [
        ("ɑ", 0.4, 0.55),
        ("a", 0.4, 0.55),
        ("ɑ", 0.0, 0.55),
        ("ɑ", 0.25, 0.55),
        ("ɑ", 0.4, 0.35),
        ("ɑ", 0.4, 0.75),
    ] {
        let (s, f) = pcm(vowel, bright, tension);
        assert!(s.iter().all(|x| x.is_finite()));
        let (band, h1h2) = spectrum(&s);
        println!(
            "{vowel} brightness={bright:.2} tension={tension:.2} H1-H2={h1h2:.2}dB bands={band:?} tract_resonances={f:?}"
        );
        assert!(band.iter().all(|v| v.is_finite()));
        assert!(f.len() >= 2 && f[0].is_finite() && f[1].is_finite());
        outcomes.push((band, h1h2, f));
    }
    // These are *relative* spectral tests, not claims about a particular person's voice.
    let main = &outcomes[0];
    let front = &outcomes[1];
    let dark = &outcomes[2];
    let lax = &outcomes[4];
    assert!(
        front.2[1] > main.2[1] + 100.0,
        "[a] should raise F2 relative to [ɑ]"
    );
    assert!(
        (front.0[2] - main.0[2]).abs() > 1.0,
        "[a] and [ɑ] need distinct resonances"
    );
    assert!(
        dark.0[4] < main.0[4] - 3.0,
        "soft closure should reduce high bands"
    );
    assert!(
        lax.1 > main.1 + 1.5,
        "relaxed phonation should increase H1-H2"
    );
    assert!(
        lax.0[4] < main.0[4] - 1.0,
        "relaxed phonation should reduce high bands"
    );
}
