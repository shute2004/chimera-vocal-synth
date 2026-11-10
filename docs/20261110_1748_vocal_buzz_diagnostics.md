# 発声の「機械的な振動音」に関する原因分析（2026-11-10）

## 観察と検証対象

利用者の試聴では、新しいネイティブマイク録音の切断感が解消した一方、以前のCHIMERA模倣候補は自然な声よりも機械的な振動音に近いと評価された。元の録音に近いピッチでも声質が一致していない。**人間の試聴を音響評価の重要な独立観測とする。**

音響診断の個人録音・個別設定・WAV・数値はAural側のGit無視対象 `out/voice-reference/` に保存する。ここでは公開できる**汎用の合成実験**だけを記録する。

## 発声源の単一要因実験

`crates/chimera/tests/voice_inverse_ablation.rs` の `glottal_pulse_shape_sensitivity` は、同じ体格・声道・圧力・ピッチを維持し、声門の開放率（open quotient）・閉鎖波形（closure sharpness）・開閉タイミング（skew）を一つずつ変更した。

| 操作 | 最終PCMのH1−H2 |
|---|---:|
| 汎用基準条件 | −10.32dB |
| open quotient = 0.45〜0.90 | −14.73〜−8.43dB |
| closure sharpness = 0〜0.9 | −10.85〜−10.12dB |
| skew = 0.35〜0.85 | −10.84〜−6.98dB |

発声源の波形によって倍音分布が変化することを確認。一方、この限定的な条件では第1倍音が第2倍音より著しく強い出力は得られていない。これを物理モデル全体の表現限界と断定することはできない。

## 声道形状の独立実験

`vowel_a_to_u_articulator_interpolation` は、CHIMERAの既存の `ArticPatch::blend` を使い、IPA `a` と `u` の**器官位置そのもの**を0〜100%の5段階で補間した。生成済みの音声をクロスフェードしたものではない。

| u側の器官姿勢 | 最終PCMのH1−H2 |
|---|---:|
| 0%（a） | −10.32dB |
| 25% | −5.85dB |
| 50% | −0.68dB |
| 75% | +6.05dB |
| 100%（u） | +23.50dB |

**声道形状だけでも低次倍音の相対強度は30dB以上変化する**。すなわち、基準録音の音色との違いを声帯発声源だけに帰することはできない。とはいえH1−H2の一致だけでは他の倍音・スペクトル包絡・声の自然さを再現できない。声道の中間姿勢、発声源、時間変動の各評価を分ける必要がある。

## 音響検査と人間の試聴の関係

Aural側で行った元録音と模倣候補のFFT比較は、ピッチが近いにもかかわらず低次倍音と非倍音成分に差が残ることを示した。粗さ・ドリフトを付与する局所実験では非倍音成分やピッチ変化の数値を調整できたが、**その数値の改善によって機械音の印象が解消されたとは確認していない**。粗さを増やすと異なる種類の不自然さが増える可能性もある。

## 検証

```sh
cd ~/Desktop/CHIMERA/chimera-vocal-synth
CARGO_INCREMENTAL=0 cargo test --offline --release -p chimera --test voice_inverse_ablation glottal_pulse_shape_sensitivity -- --ignored --nocapture
CARGO_INCREMENTAL=0 cargo test --offline --release -p chimera --test voice_inverse_ablation vowel_a_to_u_articulator_interpolation -- --ignored --nocapture
```

どちらもCHIMERA本体、既存の音声出力、母音校正テーブルを変更せず、診断テストのみ追加した。利用者の個人録音は読み込まない。
