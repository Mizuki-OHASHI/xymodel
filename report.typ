#import "@preview/physica:0.9.5": *
#import "@preview/codly:1.3.0": *

#set page(numbering: "1")
#show figure.where(kind: table): set figure.caption(position: top)

// ------------ フォント -------------
#set text(font: "Hiragino Mincho ProN")
#show heading: set text(font: "Hiragino Kaku Gothic ProN")

// ------------ コードブロック -------------
#show: codly-init.with()
#codly(
  languages: (
    gan: (name: "GaN.dat", color: color.gray),
    gaas: (name: "GaAs.dat", color: color.gray),
    alas: (name: "AlAs.dat", color: color.gray),
  ),
  zebra-fill: color.luma(248),
)

= 計算科学概論 (田浦先生担当課題)

#repeat[--]
#align(right)[
  東京大学 工学部 物理工学科4年\
  大橋 瑞輝 (学籍番号: 03-240540)\
  Email: `ohashi-mizuki0510@g.ecc.u-tokyo.ac.jp`
]
#repeat[--]

使用したPCの環境は以下の通りである。

- hardware: MacBook Air (Apple M3, arm64, 8-core CPU, 24GB RAM)
- OS: Darwin 24.5.0 (macOS 15.5)

== 問題設定

2次元XY模型のモンテカルロシミュレーションの1つであるパラレルテンパリング法 (レプリカ交換メトロポリス法) を実装し、
命令レベルの並列化、OpenMPによるスレッド並列化、SIMD命令によるデータ並列化を行い、性能を比較する#footnote[
  現在、物理工学科の授業 ($approx$ ミニ卒論) の題材として、2次元XY模型のサンプリング結果を使用する予定である。
  それに先駆けて、本課題を通して2次元XY模型のシミュレーションを行うことにした。
]。

以下、2次元XY模型およびパラレルテンパリング法の概要を、簡単に説明する。

=== 2次元XY模型

2次元XY模型は、スピン $i$ が連続的な角度 $theta_i$ を持つ2次元格子上の物理系であり、以下のハミルトニアンで定義される。
$
  H = -J sum_(angle.l i,j angle.r) cos(theta_i - theta_j)
$
J は結合定数であり、実装では $J = 1$ とする。また、$angle.l i,j angle.r$ は最近接スピンを表す。このハミルトニアンは、スピン間の相互作用を表し、スピンが同じ方向を向くほどエネルギーが低くなる。

低温においては @fig:xymodel-low のようにスピンが整列して同じ方向を向く。
一方で、高音になると @fig:xymodel-high のようにスピンの向きが乱雑になり、渦 (vortex) や反渦 (anti-vortex) が形成される。これは低音とは異なる相であり、このような相転移を KT (Kosterlitz-Thouless) 相転移と呼ぶ。

この系において重要な秩序変数 (相転移の指標) は、ヘリシティモジュラスと呼ばれる量である。
この量は以下のように定義される。

$
  Gamma
  = 1/N angle.l sum_(angle.l i,j angle.r_x) J cos(theta_i - theta_j) angle.r
  - beta/N angle.l (sum_(angle.l i,j angle.r_x) J sin(theta_i - theta_j))^2 angle.r
$

高温領域では $Gamma$ は $0$ であるが、KT 転移を境として低温側では有限の値を持つことが知られている。

本課題では、このモデルについて KT 転移が生じることを確認するために、モンテカルロシミュレーションを行い、ヘリシティモジュラスの値を計算する。なお、系のサイズは $32 times 32$ とし、周期境界条件を用いる。

#figure(caption: [2次元XY模型の例 (低温領域)])[
  #image("figs/xy_spin_field_nsamples_1000_beta_5.00_nx_32_ny_32.png", width: 75%)
]<fig:xymodel-low>

#figure(caption: [2次元XY模型の例 (高温領域)。図中の $+times$ は渦や反渦を表す。])[
  #image("figs/xy_spin_field_nsamples_1000_beta_0.80_nx_32_ny_32.png", width: 75%)
]<fig:xymodel-high>

=== メトロポリス法とパラレルテンパリング法

パラレルテンパリングは、通常のメトロポリス法を拡張したモンテカルロ手法のひとつであり、複数の温度のサンプルを並行して時間発展させ、一定のステップごとに異なる温度のサンプル間で交換を行うことで、より効率的に相空間を探索することができる。

通常のメトロポリス法ではある1つの (逆) 温度 $beta$ でシミュレーションを行う。
現在のサンプル $X_n$ があるとき、ランダムな変化 (今回の場合は特定のスピンを選んでそのスピンの角度をランダムに変化させる) を行った後のサンプル $X'$ を作る。
$X_n$ と $X'$ のエネルギーの差を $Delta E = E(X') - E(X_n)$ とすると、$min(1, exp(- beta Delta E))$ の確率で $X'$ を受け入れて $X_{n+1} <- X'$ とする。
こうすることで $X_n$ と $X_(n+1)$ の間に詳細釣り合いの条件が成り立つので、十分長い時間このステップを続けると、サンプルは熱平衡状態 (カノニカル分布) に収束する。

パラレルテンパリングでは、複数の温度 $beta_1, beta_2, ..., beta_N$ のサンプルを用意し、それぞれのサンプルに対してメトロポリス法を適用し、発展の途中の一定の間隔で、異なる温度のサンプル間で交換を行う。
この際にはふたつの温度のサンプルに関して温度差を $Delta beta$、エネルギー差を $Delta E$ とすると、交換を行う確率は $min(1, exp(Delta beta Delta E))$ である。
こうすることで、異なる温度のサンプル間でも詳細釣り合いの条件が成り立つので、全てのサンプルは熱平衡状態に収束する。

== 実装

ソースコードは #link("git@github.com:Mizuki-OHASHI/xymodel.git")[`GitHub:Mizuki-OHASHI/xymodel`] で公開しているが、本レポートの末尾にもソースコードの一部を掲載する。

ベースとなるコードは、上の XY 模型の説明をプロンプトに付した上で Gemini に生成させた C 言語のコードである (Appendix 1)。上で説明したパラレルテンパリング法が実装されており、各温度でのヘリシティモジュラスのサンプリング平均や実行時間 (wall time) を出力する。

ここでは、並列化の理解のために、処理の流れを簡単な疑似コードで示す。

```algorithm
時間の計測開始
各温度のサンプリングを初期化する
物理量_total = 0.0
i = 0
while (サンプリングした数 < サンプリング数) {
  for (各温度のサンプリング) { ・・・(1)
    全てのサイトのスピンについてメトロポリス法を適用する ・・・(2)
  }
  for (複数回) { ・・・(3)
    隣り合った温度のサンプリング間でレプリカ交換を行う ・・・(4)
  }
  if (初期緩和後のサンプリング間隔に達したら) { ・・・(5)
    各温度のサンプリングから物理量を測定する ・・・(6)
    物理量_total += 測定した物理量
  }
  i ++
}
出力ファイルに物理量の平均値 (物理量_total / サンプリング数) を書き込む
時間の計測終了
```

実行結果を @fig:helicity-modulus-original-noopt に示す。
実行時間は *$135$ 秒*であった。
プロットを見ると、$beta = 0.9$ の付近でヘリシティモジュラスが $0$ から有限の値に変化していることがわかり、KT 転移が確認できる。

なお、無限大の格子サイズにおける転移 (逆) 温度は $1.1$ 程度であることが知られているが、格子サイズが有限である #footnote[2 次元 XY モデルは相関長が無限大に発散することから、サイズが有限であることの影響を大きく受ける。] ことから、転移温度が高温側にずれていると考えられる。

#figure(caption: [ヘリシティモジュラスのサンプリング平均 (オリジナルのコード)\
  (コンパイル: `clang original.c -o original -lm`)])[
  #image("figs/original_noopt.dat.png", width: 75%)
]<fig:helicity-modulus-original-noopt>

以下、並列化などの工夫を行なった場合の結果が正しいことを確認する際には、
(乱数を含む計算のため) 完全に同じ結果が得られないため、このプロットと同様のデータが得られることを確認する。

=== 命令レベルの並列化

ソースコード: 変更なし

まずは、コンパイル時に最適化オプションを指定して、命令レベルの並列化を行う。
`clang -O3 original.c -o original_noopt -lm` としてコンパイルしたコードを実行すると、得られる結果は @fig:helicity-modulus-original-noopt と全く同じで、実行時間は *$97.2$ 秒* まで短縮された (最初と比較して 0.72 倍) 。

=== OpenMPによるスレッド並列化

ソースコード: `openmp.c` (GitHub)

続いて、OpenMPを用いてスレッド並列化を行う。
具体的には上の疑似コードにおいて以下の 2 箇所のループを並列化する。

- `(2)` のメトロポリス法の処理が重いので、各温度のサンプリングを並列化する (ループ `(1)` の並列化)
- `(6)` の物理量の計算は各々独立に実行可能なので、各温度のサンプリングごとに並列化する (ループ `(5)` の並列化)

なお、並列化を行わなかった箇所があるが、それは以下の理由による。

- `(2)` のメトロポリス法の処理の内部については、並列化の余地はあるものの、素朴にスピンごとの更新を並列化すると、隣り合ったスピンの更新が競合してしまう恐れがあるので、今は並列化しない。
- `(4)` のレプリカ交換については、アドレスの入れ替えをしているだけなので各処理はそこまで重くなく、しかも並列化することで隣り合ったレプリカ同士の競合が発生する可能性があるため、並列化しない。

なお、標準的に用意されている乱数生成器は、グルーバルに状態を保つことから並列化すると競合が発生する恐れがあるため、
スレッドごとに独立な乱数生成器を用意して、各スレッドで独立に乱数を生成するようにした。
そのためにオープンソースで公開されている #link("https://github.com/imneme/pcg-c/tree/master")[PCG Random Number Generation, C Edition] を用いた。

以上を並列化を行なって得られた結果を @fig:helicity-modulus-openmp に示す。
実行時間は *$25.3$ 秒*であった (最初と比較して 0.19 倍)。
プロットから、正しくシミュレーションが行われていると考えられる。

#figure(caption: [ヘリシティモジュラスのサンプリング平均 (OpenMPによる並列化)\
  (コンパイル: `clang -Xpreprocessor -fopenmp -Rpass=loop-vectorize -march=native -ffast-math -lomp -I"$(shell brew --prefix libomp)/include" -L"$(shell brew --prefix libomp)/lib" -O3 simd.c -o simd -lm`)])[
  #image("figs/openmp.dat.png", width: 75%)
]<fig:helicity-modulus-openmp>

=== SIMD命令によるデータ並列化

ソースコード: `simd.c` (GitHub または Appendix 2)

次に上の並列化に加えて、SIMD命令を用いてデータ並列化を行う。

`(2)` メトロポリス法の処理は、素朴に並列化すると、隣同士のスピンを同時に更新して競合してしまう
(隣のスピンに依存して更新が行われるので、同時に更新されると不適当な値が計算される可能性がある)。
そこで、スピンの更新順序を工夫することで競合を回避する。
具体的な工夫は @fig:spin-update-order の通りである。

オリジナルのコードでは、左上から順番にひとつずつスピンを更新していた。
そのため、そのままの順番で並列化すると、隣り合ったスピンの更新が同期してしまう。
そこで、市松模様上にスピンを更新するように順序を変更した。スピンの更新が依存するのは、上下左右の隣接スピンのみであるため、市松模様上のスピンは互いに独立に更新できる。

#figure(caption: [スピンの更新順序の工夫])[
  #image("figs/ichimatsu.jpeg", width: 75%)
]<fig:spin-update-order>

順序を変更した上で、

+ 1 個飛ばしのスピンを連続するメモリ領域にコピーして (gather load)
+ IMD 命令によって一括で更新する。
+ 最後にメトロポリス法の判定で受け入れられたスピンを更新する (scatter store)

このようにすることで、スピンの更新を SIMD 命令で並列化することができる。
該当箇所のみ、ソースコードを示すと以下のようになる。

#text(size: 8pt)[
  ```c
  // メトロポリス法の SIMD 化
  void metropolis_sweep(Replica *rep, pcg32_random_t *rng)
  {
    const int N_SITES = rep->nx * rep->ny;
    const float DELTA = 1.0f; // スピン角度の更新幅

    int odd, site_idx, idx, ix, iy, r; // loop variables

    // block variables for SIMD
    float spin[NL], new_spin[NL];                                         // spin block
    float neighborsR[NL], neighborsL[NL], neighborsD[NL], neighborsU[NL]; // neighbors blocks
    float delta_e[NL];                                                    // energy change and acceptance probability
    float rand_val_acc[NL], rand_val_delta[NL];                           // random values for acceptance and delta
    int accept[NL];                                                       // acceptance flags

    // update each site in blocks of NL
    for (odd = 0; odd < 2; odd++)
    {
      // odd = 0 : even sites
      // odd = 1 : odd sites
      for (site_idx = odd; site_idx < N_SITES; site_idx += NL * 2)
      {
        // gather load
        for (r = 0; r < NL; ++r)
        {
          idx = site_idx + r * 2;
          ix = idx % rep->nx;
          iy = idx / rep->nx;

          // load spin
          spin[r] = rep->spin[idx];

          // load neighbors
          neighborsR[r] = rep->spin[iy * rep->nx + (ix + 1) % rep->nx];             // Right
          neighborsL[r] = rep->spin[iy * rep->nx + (ix - 1 + rep->nx) % rep->nx];   // Left
          neighborsD[r] = rep->spin[((iy + 1) % rep->ny) * rep->nx + ix];           // Down
          neighborsU[r] = rep->spin[((iy - 1 + rep->ny) % rep->ny) * rep->nx + ix]; // Up

          rand_val_delta[r] = (float)pcg32_random_r(rng) / (float)UINT32_MAX;
          rand_val_acc[r] = (float)pcg32_random_r(rng) / (float)UINT32_MAX;
        }

  #pragma omp simd
        for (int r = 0; r < NL; ++r)
        {
          new_spin[r] = fmodf(spin[r] + (rand_val_delta[r] - 0.5f) * DELTA + 2.0f * M_PI, 2.0f * M_PI);

          // calculate energy change
          delta_e[r] = 0.0f;
          // calculate energy change with neighbors
          delta_e[r] += -J * (cosf(new_spin[r] - neighborsR[r]) - cosf(spin[r] - neighborsR[r]));
          delta_e[r] += -J * (cosf(new_spin[r] - neighborsL[r]) - cosf(spin[r] - neighborsL[r]));
          delta_e[r] += -J * (cosf(new_spin[r] - neighborsD[r]) - cosf(spin[r] - neighborsD[r]));
          delta_e[r] += -J * (cosf(new_spin[r] - neighborsU[r]) - cosf(spin[r] - neighborsU[r]));

          accept[r] = (delta_e[r] < 0.0f || rand_val_acc[r] < expf(-rep->beta * delta_e[r]));
        }

        // scatter store
        for (int r = 0; r < NL; ++r)
        {
          int idx = site_idx + r * 2;
          if (accept[r])
          {
            // accept the new spin
            rep->spin[idx] = new_spin[r];
            rep->energy += delta_e[r];
          }
        }
      }
    }
  }
  ```
]

このSIMD化を行った上で、`clang -Xpreprocessor -fopenmp -Rpass=loop-vectorize -march=native -ffast-math -lomp -I"$(shell brew --prefix libomp)/include" -L"$(shell brew --prefix libomp)/lib" -O3 simd.c -o simd -lm` としてコンパイルした。
まず、コンパイル時に出た以下のメッセージから、正常に SIMD 化が行われたことがわかる。

```
simd.c:143:1: remark: vectorized loop (vectorization width: 4, interleaved count: 1) [-Rpass=loop-vectorize]
  143 | #pragma omp simd
      | ^
simd.c:318:1: remark: vectorized loop (vectorization width: 2, interleaved count: 4) [-Rpass=loop-vectorize]
  318 | #pragma omp simd reduction(+ : E_x_sweep, J_x_sweep, E_y_sweep, J_y_sweep)
      | ^
```

得られた結果を @fig:helicity-modulus-simd に示す。
結果はオリジナルのソースコードと同様に KT 転移が確認でき、問題なくシミュレーションが行われていることがわかる。

#figure(caption: [ヘリシティモジュラスのサンプリング平均 (SIMD命令によるデータ並列化)\
  (コンパイル: `clang -Xpreprocessor -fopenmp -Rpass=loop-vectorize -march=native -ffast-math -lomp -I"$(shell brew --prefix libomp)/include" -L"$(shell brew --prefix libomp)/lib" -O3 simd.c -o simd -lm`)])[
  #image("figs/simd.dat.png", width: 75%)
]<fig:helicity-modulus-simd>

実行時間は *$30.0$ 秒* であった (最初と比較して 0.22 倍)。
SIMD化を行ったことでむしろ実行時間が増加してしまい、並列化の効果よりオーバーヘッドの方が大きくなってしまった。
これには、いくつかの要因が考えられるので、以下に挙げる。

- SIMD命令による一括計算をする工程があまり多くなく (エネルギーの差分を計算するだけ) 、並列化の効果が薄いこと。
- SIMD命令を用いるために連続したメモリ領域にスピンのデータを複製した (gather load、上のコードの 23 行目以降のブロック) が、このオーバーヘッドが大きいこと。
- メトロポリス法の性質上、最後に受容確率に応じてスピンを更新する・しないの `IF` 判定が必須である。実装では scatter store の際にこの判定を行うことで、メインの計算部分を一括化しているが、この判定が並列化の阻害要因となって、十分に並列化できなかった。

これらについては、より詳細に原因を調査し、改善する余地があると考えられる。
さらに高速化する提案を以下に挙げる。

==== Gather load のオーバーヘッドとスピンの保持の仕方

ひとつ飛ばし (市松模様) のスピンを連続したメモリ領域にコピーする (gather load) のオーバーヘッドが大きいことが原因の一つであると考えられる。そのため、そもそものスピン情報の格納方法を工夫することで、gather load のオーバーヘッドを減らすことができる。
つまり、ひとつ飛ばしのスピンが元々連続したメモリ領域に格納されれば良い。

さらに、周期境界条件を考慮するにあたって、現在は例えば $x$ 方向の隣のスピンを取得する際に
```
spin[iy * nx + (ix + 1) % nx]
```
のように、`%` 演算子を用いて周期境界条件を適用しているが、これもオーバーヘッドとなる。
なぜなら、アドレス計算のたびに `%` 演算を行う負荷があるだけでなく、メモリが不連続になるからである。
そこで、境界より一回りだけ大きいサイズのメモリ領域を確保して、境界の値が反対の端と同じになるように管理する。

そのようにすることで、境界のスピンを取得する際に `%` 演算を行う必要がなり、境界付近も含めてスピンのデータを連続したメモリ領域に格納することができる (@fig:spin-layout)。

#figure(caption: [スピンの格納方法の工夫\
  #align(left)[#text(size: 10pt)[
      境界部分のセルを一回り大きく拡張する。
      この際周期境界条件に注意する (図においては、$star$ マーク同士のセルは等価なセルとして扱い、同じ値を持つように管理する)。このように拡張することで、境界のスピンの更新を行う際にも、境界から離れたスピンと同様に、`spin[iy * nx + (ix + 1)]` のように周期境界条件を考慮せずに取得できる
      (図においては セル1 の更新をするときは白抜きされた赤丸を参照すればよく、
      塗りつぶされた赤丸のセルを参照しなくて良くなる。
      これによって `%` 演算が不要になるとともに、 セル 1,2, ... を並列に計算するにあたって、対象のセルの隣のセル同士もメモリ領域において必ず連続するようになり、効率的な並列計算が可能になると考えられる)。
    ]]])[
  #image("figs/memory_layout.jpeg", width: 50%)
]<fig:spin-layout>

==== 条件分岐の回避 (組み込み関数の利用)

メトロポリス法の受容確率の判定において、`IF` 文を用いてスピンを更新する・しないを決定しているが、この条件分岐が並列化の阻害要因となっている。
組み込み関数 mask を用いることで、この条件分岐を回避することができそうである。
調べたところによると、ハードウェア環境の依存なども大きいようであり、今回は断念した。

== Appendix: ソースコード

=== Appendix 1. 並列化等の工夫をしていないオリジナルのコード

ファイル名: `original.c`

#text(size: 8pt)[
  ```c
  // オリジナルの XY モデルのパラレルテンパリング (レプリカ交換) シミュレーション
  // Gemini によって生成されたコードをベースとして, 部分的に修正を加えた (目安としてコメントが英語の部分は修正した箇所)

  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <time.h>
  #include <sys/time.h>

  #define NX 32                       // 格子のXサイズ
  #define NY 32                       // 格子のYサイズ
  #define N_REPLICAS 16               // レプリカの数
  #define N_SAMPLING 1000             // サンプリング数
  #define SAMPLE_INTERVAL_SWEEPS 100  // サンプリング間隔
  #define THERMALIZATION_SWEEPS 10000 // 初期緩和のスイープ数

  // 相互作用の強さ
  #define J (1.0f)

  /**
   * @brief レプリカの状態を保持する構造体
   */
  typedef struct
  {
    float *spin;  // スピンの角度θを格納する配列 (サイズ: N_SITES)
    float beta;   // このレプリカの逆温度 β = 1/T
    float energy; // このレプリカの現在のエネルギー
    int nx;       // 格子のX方向のサイズ
    int ny;       // 格子のY方向のサイズ
  } Replica;

  /**
   * @brief 系全体のエネルギーを計算する
   * @param rep 計算対象のレプリカへのポインタ
   * @return 計算された全エネルギー
   * @note 周期境界条件を適用します
   */
  float calculate_total_energy(Replica *rep)
  {
    float total_energy = 0.0f;
    const int N_SITES = rep->nx * rep->ny;

    for (int i = 0; i < N_SITES; ++i)
    {
      int ix = i % rep->nx;
      int iy = i / rep->nx;

      // 右の隣接スピンとの相互作用
      int right_neighbor_idx = iy * rep->nx + (ix + 1) % rep->nx;
      total_energy -= J * cosf(rep->spin[i] - rep->spin[right_neighbor_idx]);

      // 下の隣接スピンとの相互作用
      int down_neighbor_idx = ((iy + 1) % rep->ny) * rep->nx + ix;
      total_energy -= J * cosf(rep->spin[i] - rep->spin[down_neighbor_idx]);
    }
    return total_energy;
  }

  /**
   * @brief パラレルテンパレリング用のレプリカ群を初期化する
   * @param replicas レプリカの配列へのポインタ
   * @param n_replicas レプリカの数
   * @param nx 格子のX方向のサイズ
   * @param ny 格子のY方向のサイズ
   * @param betas 各レプリカに設定する逆温度の配列
   * @note 乱数シードは事前に設定しておく必要があります
   */
  void initialize_replicas(Replica *replicas, int n_replicas, int nx, int ny, const float *betas)
  {
    const int N_SITES = nx * ny;
    for (int i = 0; i < n_replicas; ++i)
    {
      replicas[i].spin = (float *)malloc(sizeof(float) * N_SITES);
      if (replicas[i].spin == NULL)
      {
        fprintf(stderr, "Error: Failed to allocate memory for spins.\n");
        exit(EXIT_FAILURE);
      }

      replicas[i].beta = betas[i];
      replicas[i].nx = nx;
      replicas[i].ny = ny;

      // スピンをランダムな角度 [0, 2π) で初期化
      for (int j = 0; j < N_SITES; ++j)
      {
        replicas[i].spin[j] = ((float)rand() / (float)RAND_MAX) * 2.0f * M_PI;
      }

      // 初期エネルギーを計算
      replicas[i].energy = calculate_total_energy(&replicas[i]);
    }
  }

  /**
   * @brief 1つのレプリカに対して1モンテカルロスイープを実行する (メトロポリス法)
   * @param rep 更新するレプリカへのポインタ
   * @note 全てのスピンサイトを0からN_SITES-1まで順番に1回ずつ更新します。
   */
  void metropolis_sweep(Replica *rep)
  {
    const int N_SITES = rep->nx * rep->ny;
    const float DELTA = 1.0f; // スピン角度の更新幅

    // 全てのサイトを0からN_SITES-1まで順番にループ
    for (int site_idx = 0; site_idx < N_SITES; ++site_idx)
    {
      int ix = site_idx % rep->nx;
      int iy = site_idx / rep->nx;

      float old_spin = rep->spin[site_idx];
      // 乱数生成器はスピンの新しい角度を試すために使用
      float new_spin = fmodf(old_spin + ((float)rand() / (float)RAND_MAX - 0.5f) * DELTA, 2.0f * M_PI);
      if (new_spin < 0.0f)
      {
        new_spin += 2.0f * M_PI;
      }

      // エネルギー変化量を計算
      float delta_e = 0.0f;
      // 4つの隣接サイトをループ
      int neighbors[4];
      neighbors[0] = iy * rep->nx + (ix + 1) % rep->nx;             // Right
      neighbors[1] = iy * rep->nx + (ix - 1 + rep->nx) % rep->nx;   // Left
      neighbors[2] = ((iy + 1) % rep->ny) * rep->nx + ix;           // Down
      neighbors[3] = ((iy - 1 + rep->ny) % rep->ny) * rep->nx + ix; // Up

      for (int j = 0; j < 4; ++j)
      {
        int neighbor_idx = neighbors[j];
        // new_spinと隣接スピンとの相互作用エネルギーと、
        // old_spinと隣接スピンとの相互作用エネルギーの差を計算
        delta_e += -J * (cosf(new_spin - rep->spin[neighbor_idx]) - cosf(old_spin - rep->spin[neighbor_idx]));
      }

      // メトロポリス判定
      if (delta_e < 0.0f || ((float)rand() / (float)RAND_MAX) < expf(-rep->beta * delta_e))
      {
        rep->spin[site_idx] = new_spin;
        rep->energy += delta_e;
      }
    }
  }

  /**
   * @brief レプリカ交換を実行する
   * @param replicas レプリカの配列
   * @param n_replicas レプリカの数
   */
  void replica_exchange(Replica *replicas, int n_replicas)
  {
    // 隣接するレプリカのペア (i, i+1) をランダムに選択
    int i = (rand() % (n_replicas - 1));

    float delta_beta = replicas[i].beta - replicas[i + 1].beta;
    float delta_energy = replicas[i].energy - replicas[i + 1].energy;

    // 交換確率を計算: min(1, exp(Δβ * ΔE))
    float acceptance_prob = expf(delta_beta * delta_energy);

    if (rand() < (float)RAND_MAX * acceptance_prob)
    {
      // スピン配列とエネルギーを交換
      float *temp_spin = replicas[i].spin;
      replicas[i].spin = replicas[i + 1].spin;
      replicas[i + 1].spin = temp_spin;

      float temp_energy = replicas[i].energy;
      replicas[i].energy = replicas[i + 1].energy;
      replicas[i + 1].energy = temp_energy;
    }
  }

  /**
   * @brief 確保したメモリを解放する
   * @param replicas レプリカの配列
   * @param n_replicas レプリカの数
   */
  void free_replicas(Replica *replicas, int n_replicas)
  {
    for (int i = 0; i < n_replicas; ++i)
    {
      if (replicas[i].spin != NULL)
      {
        free(replicas[i].spin);
        replicas[i].spin = NULL;
      }
    }
  }

  int main(int argc, char *argv[])
  {
    // 1. パラメータ設定
    const int N_SITES = NX * NY;
    Replica replicas[N_REPLICAS];
    float betas[N_REPLICAS];
    int i, sweep;

    float min_beta = 0.4f;
    float max_beta = 1.9f;

    struct timeval start_time, end_time;
    double wall_time;

    FILE *fp; // file pointer for output

    // start time measurement
    gettimeofday(&start_time, NULL);

    if (argc != 2)
    {
      fprintf(stderr, "Usage: %s <output_file>\n", argv[0]);
      return EXIT_FAILURE;
    }

    // 出力ファイルを開く
    fp = fopen(argv[1], "w");
    if (fp == NULL)
    {
      fprintf(stderr, "Error: Could not open output file %s\n", argv[1]);
      return EXIT_FAILURE;
    }

    fprintf(fp, "<version\noriginal\nversion>\n");

    fprintf(fp, "<input_parameters\n");
    fprintf(fp, "N_SITES=%d N_REPLICAS=%d NX=%d NY=%d N_SAMPLING=%d SAMPLE_INTERVAL_SWEEPS=%d THERMALIZATION_SWEEPS=%d\n",
            N_SITES, N_REPLICAS, NX, NY, N_SAMPLING, SAMPLE_INTERVAL_SWEEPS, THERMALIZATION_SWEEPS);
    fprintf(fp, "min_beta=%.6f max_beta=%.6f\n", min_beta, max_beta);
    fprintf(fp, "input_parameters>\n");

    for (i = 0; i < N_REPLICAS; ++i)
    {
      betas[i] = min_beta + (max_beta - min_beta) * i / (N_REPLICAS - 1);
    }

    // 2. 初期化
    srand(48);
    initialize_replicas(replicas, N_REPLICAS, NX, NY, betas);

    // 全レプリカの物理量を格納する配列
    double total_Ex[N_REPLICAS] = {0.0};
    double total_Jx_squared[N_REPLICAS] = {0.0};
    double total_Ey[N_REPLICAS] = {0.0};
    double total_Jy_squared[N_REPLICAS] = {0.0};
    long measurement_count = 0;

    // 3. シミュレーションループ
    printf("<simulation_progress\n");
    sweep = 0;
    while (measurement_count < N_SAMPLING)
    {
      // 各レプリカでメトロポリス更新
      for (i = 0; i < N_REPLICAS; ++i)
      {
        metropolis_sweep(&replicas[i]);
      }

      // レプリカ交換
      for (i = 0; i < N_REPLICAS; ++i)
      { // 交換頻度を上げるためにループ
        replica_exchange(replicas, N_REPLICAS);
      }

      // 物理量の測定(初期緩和後)
      if (sweep > THERMALIZATION_SWEEPS && (sweep - THERMALIZATION_SWEEPS) % SAMPLE_INTERVAL_SWEEPS == 0)
      {
        for (i = 0; i < N_REPLICAS; ++i)
        {
          float E_x_sweep = 0.0f;
          float J_x_sweep = 0.0f;
          float E_y_sweep = 0.0f;
          float J_y_sweep = 0.0f;

          for (int j = 0; j < N_SITES; ++j)
          {
            int ix = j % NX;
            int iy = j / NX;

            // x方向 (右の隣人)
            int right_neighbor_idx = iy * NX + (ix + 1) % NX;
            float delta_theta_x = replicas[i].spin[j] - replicas[i].spin[right_neighbor_idx];
            E_x_sweep += J * cosf(delta_theta_x);
            J_x_sweep += J * sinf(delta_theta_x);

            // y方向 (下の隣人)
            int down_neighbor_idx = ((iy + 1) % NY) * NX + ix;
            float delta_theta_y = replicas[i].spin[j] - replicas[i].spin[down_neighbor_idx];
            E_y_sweep += J * cosf(delta_theta_y);
            J_y_sweep += J * sinf(delta_theta_y);
          }
          total_Ex[i] += E_x_sweep;
          total_Jx_squared[i] += (double)J_x_sweep * J_x_sweep;
          total_Ey[i] += E_y_sweep;
          total_Jy_squared[i] += (double)J_y_sweep * J_y_sweep;
        }
        measurement_count++;
        if (measurement_count % (N_SAMPLING / 10) == 0)
        {
          printf("sweep=%d measurement_count=%ld total_Ex=%f total_Jx_squared=%f total_Ey=%f total_Jy_squared=%f\n",
                 sweep, measurement_count, total_Ex[0], total_Jx_squared[0], total_Ey[0], total_Jy_squared[0]);
        }
      }

      sweep++;
    }
    printf("simulation_progress>\n");

    // 最終結果の計算と表示
    fprintf(fp, "<helicity_modulus_results\n");
    for (i = 0; i < N_REPLICAS; ++i)
    {
      if (measurement_count > 0)
      {
        double avg_Ex = total_Ex[i] / measurement_count;
        double avg_Jx_squared = total_Jx_squared[i] / measurement_count;
        double upsilon_x = (avg_Ex / N_SITES) - (replicas[i].beta / N_SITES) * avg_Jx_squared;

        double avg_Ey = total_Ey[i] / measurement_count;
        double avg_Jy_squared = total_Jy_squared[i] / measurement_count;
        double upsilon_y = (avg_Ey / N_SITES) - (replicas[i].beta / N_SITES) * avg_Jy_squared;

        double helicity_modulus = (upsilon_x + upsilon_y) / 2.0;
        double energy_per_site = replicas[i].energy / (double)N_SITES;
        fprintf(fp, "replica=%d beta=%.6f energy_per_site=%.6f helicity_modulus=%.6f\n",
                i, replicas[i].beta, energy_per_site, helicity_modulus);
      }
    }
    fprintf(fp, "helicity_modulus_results>\n");

    // 4. メモリ解放
    free_replicas(replicas, N_REPLICAS);

    // end time measurement
    gettimeofday(&end_time, NULL);
    wall_time = (end_time.tv_sec - start_time.tv_sec) +
                (end_time.tv_usec - start_time.tv_usec) / 1000000.0;
    fprintf(fp, "<execution_time\n%.4f\nexecution_time>\n", wall_time);

    fclose(fp);

    return 0;
  }
  ```
]

=== Appendix 2. SIMD命令によるデータ並列化を施したコード

ファイル名: `simd.c`

#text(size: 8pt)[
  ```c
  // CPUでの性能向上 --- マルチコア向上 (OpenMP) + SIMD 化

  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <time.h>
  #include <sys/time.h>
  #include <omp.h>
  #include "pcg_variants.h"

  #define NX 32                       // 格子のXサイズ
  #define NY 32                       // 格子のYサイズ
  #define NL 16                       // SIMD化のブロックサイズ
  #define N_REPLICAS 16               // レプリカの数
  #define N_SAMPLING 1000             // サンプリング数
  #define SAMPLE_INTERVAL_SWEEPS 100  // サンプリング間隔
  #define THERMALIZATION_SWEEPS 10000 // 初期緩和のスイープ数

  // 相互作用の強さ
  #define J (1.0f)

  /**
   * @brief レプリカの状態を保持する構造体
   */
  typedef struct
  {
    float *spin;  // スピンの角度θを格納する配列 (サイズ: N_SITES)
    float beta;   // このレプリカの逆温度 β = 1/T
    float energy; // このレプリカの現在のエネルギー
    int nx;       // 格子のX方向のサイズ
    int ny;       // 格子のY方向のサイズ
  } Replica;

  /**
   * @brief 系全体のエネルギーを計算する
   * @param rep 計算対象のレプリカへのポインタ
   * @return 計算された全エネルギー
   * @note 周期境界条件を適用します
   */
  float calculate_total_energy(Replica *rep)
  {
    float total_energy = 0.0f;
    const int N_SITES = rep->nx * rep->ny;

    for (int i = 0; i < N_SITES; ++i)
    {
      int ix = i % rep->nx;
      int iy = i / rep->nx;

      // 右の隣接スピンとの相互作用
      int right_neighbor_idx = iy * rep->nx + (ix + 1) % rep->nx;
      total_energy -= J * cosf(rep->spin[i] - rep->spin[right_neighbor_idx]);

      // 下の隣接スピンとの相互作用
      int down_neighbor_idx = ((iy + 1) % rep->ny) * rep->nx + ix;
      total_energy -= J * cosf(rep->spin[i] - rep->spin[down_neighbor_idx]);
    }
    return total_energy;
  }

  /**
   * @brief パラレルテンパレリング用のレプリカ群を初期化する
   * @param replicas レプリカの配列へのポインタ
   * @param n_replicas レプリカの数
   * @param nx 格子のX方向のサイズ
   * @param ny 格子のY方向のサイズ
   * @param betas 各レプリカに設定する逆温度の配列
   * @note 乱数シードは事前に設定しておく必要があります
   */
  void initialize_replicas(Replica *replicas, int n_replicas, int nx, int ny, const float *betas)
  {
    const int N_SITES = nx * ny;
    for (int i = 0; i < n_replicas; ++i)
    {
      replicas[i].spin = (float *)malloc(sizeof(float) * N_SITES);
      if (replicas[i].spin == NULL)
      {
        fprintf(stderr, "Error: Failed to allocate memory for spins.\n");
        exit(EXIT_FAILURE);
      }

      replicas[i].beta = betas[i];
      replicas[i].nx = nx;
      replicas[i].ny = ny;

      // スピンをランダムな角度 [0, 2π) で初期化
      for (int j = 0; j < N_SITES; ++j)
      {
        replicas[i].spin[j] = ((float)rand() / (float)RAND_MAX) * 2.0f * M_PI;
      }

      // 初期エネルギーを計算
      replicas[i].energy = calculate_total_energy(&replicas[i]);
    }
  }

  /**
   * @brief 1つのレプリカに対して1モンテカルロスイープを実行する (メトロポリス法)
   * @param rep 更新するレプリカへのポインタ
   * @note 全てのスピンサイトを更新します
   */
  void metropolis_sweep(Replica *rep, pcg32_random_t *rng)
  {
    const int N_SITES = rep->nx * rep->ny;
    const float DELTA = 1.0f; // スピン角度の更新幅

    int odd, site_idx, idx, ix, iy, r; // loop variables

    // block variables for SIMD
    float spin[NL], new_spin[NL];                                         // spin block
    float neighborsR[NL], neighborsL[NL], neighborsD[NL], neighborsU[NL]; // neighbors blocks
    float delta_e[NL];                                                    // energy change and acceptance probability
    float rand_val_acc[NL], rand_val_delta[NL];                           // random values for acceptance and delta
    int accept[NL];                                                       // acceptance flags

    // update each site in blocks of NL
    for (odd = 0; odd < 2; odd++)
    {
      // odd = 0 : even sites
      // odd = 1 : odd sites
      for (site_idx = odd; site_idx < N_SITES; site_idx += NL * 2)
      {
        // gather load
        for (r = 0; r < NL; ++r)
        {
          idx = site_idx + r * 2;
          ix = idx % rep->nx;
          iy = idx / rep->nx;

          // load spin
          spin[r] = rep->spin[idx];

          // load neighbors
          neighborsR[r] = rep->spin[iy * rep->nx + (ix + 1) % rep->nx];             // Right
          neighborsL[r] = rep->spin[iy * rep->nx + (ix - 1 + rep->nx) % rep->nx];   // Left
          neighborsD[r] = rep->spin[((iy + 1) % rep->ny) * rep->nx + ix];           // Down
          neighborsU[r] = rep->spin[((iy - 1 + rep->ny) % rep->ny) * rep->nx + ix]; // Up

          rand_val_delta[r] = (float)pcg32_random_r(rng) / (float)UINT32_MAX;
          rand_val_acc[r] = (float)pcg32_random_r(rng) / (float)UINT32_MAX;
        }

        // scatter store
  #pragma omp simd
        for (int r = 0; r < NL; ++r)
        {
          new_spin[r] = fmodf(spin[r] + (rand_val_delta[r] - 0.5f) * DELTA + 2.0f * M_PI, 2.0f * M_PI);

          // calculate energy change
          delta_e[r] = 0.0f;
          // calculate energy change with neighbors
          delta_e[r] += -J * (cosf(new_spin[r] - neighborsR[r]) - cosf(spin[r] - neighborsR[r]));
          delta_e[r] += -J * (cosf(new_spin[r] - neighborsL[r]) - cosf(spin[r] - neighborsL[r]));
          delta_e[r] += -J * (cosf(new_spin[r] - neighborsD[r]) - cosf(spin[r] - neighborsD[r]));
          delta_e[r] += -J * (cosf(new_spin[r] - neighborsU[r]) - cosf(spin[r] - neighborsU[r]));

          accept[r] = (delta_e[r] < 0.0f || rand_val_acc[r] < expf(-rep->beta * delta_e[r]));
        }

        // store results back to replica
        for (int r = 0; r < NL; ++r)
        {
          int idx = site_idx + r * 2;
          if (accept[r])
          {
            // accept the new spin
            rep->spin[idx] = new_spin[r];
            rep->energy += delta_e[r];
          }
        }
      }
    }
  }

  /**
   * @brief レプリカ交換を実行する
   * @param replicas レプリカの配列
   * @param n_replicas レプリカの数
   */
  void replica_exchange(Replica *replicas, int n_replicas)
  {
    // 隣接するレプリカのペア (i, i+1) をランダムに選択
    int i = (rand() % (n_replicas - 1));

    float delta_beta = replicas[i].beta - replicas[i + 1].beta;
    float delta_energy = replicas[i].energy - replicas[i + 1].energy;

    // 交換確率を計算: min(1, exp(Δβ * ΔE))
    float acceptance_prob = expf(delta_beta * delta_energy);

    if (rand() < (float)RAND_MAX * acceptance_prob)
    {
      // スピン配列とエネルギーを交換
      float *temp_spin = replicas[i].spin;
      replicas[i].spin = replicas[i + 1].spin;
      replicas[i + 1].spin = temp_spin;

      float temp_energy = replicas[i].energy;
      replicas[i].energy = replicas[i + 1].energy;
      replicas[i + 1].energy = temp_energy;
    }
  }

  /**
   * @brief 確保したメモリを解放する
   * @param replicas レプリカの配列
   * @param n_replicas レプリカの数
   */
  void free_replicas(Replica *replicas, int n_replicas)
  {
    for (int i = 0; i < n_replicas; ++i)
    {
      if (replicas[i].spin != NULL)
      {
        free(replicas[i].spin);
        replicas[i].spin = NULL;
      }
    }
  }

  int main(int argc, char *argv[])
  {
    // 1. パラメータ設定
    const int N_SITES = NX * NY;
    Replica replicas[N_REPLICAS];
    float betas[N_REPLICAS];
    int i, sweep;

    float min_beta = 0.4f;
    float max_beta = 1.9f;

    struct timeval start_time, end_time;
    double wall_time;

    FILE *fp;

    // start time measurement
    gettimeofday(&start_time, NULL);

    if (N_SITES % (NL * 2) != 0)
    {
      fprintf(stderr, "Error: N_SITES must be a multiple of %d.\n", NL * 2);
      return EXIT_FAILURE;
    }

    if (argc != 2)
    {
      fprintf(stderr, "Usage: %s <output_file>\n", argv[0]);
      return EXIT_FAILURE;
    }

    // 出力ファイルを開く
    fp = fopen(argv[1], "w");
    if (fp == NULL)
    {
      fprintf(stderr, "Error: Could not open output file %s\n", argv[1]);
      return EXIT_FAILURE;
    }

    fprintf(fp, "<version\nsimd\nversion>\n");

    fprintf(fp, "<input_parameters\n");
    fprintf(fp, "N_SITES=%d N_REPLICAS=%d NX=%d NY=%d N_SAMPLING=%d SAMPLE_INTERVAL_SWEEPS=%d THERMALIZATION_SWEEPS=%d\n",
            N_SITES, N_REPLICAS, NX, NY, N_SAMPLING, SAMPLE_INTERVAL_SWEEPS, THERMALIZATION_SWEEPS);
    fprintf(fp, "min_beta=%.6f max_beta=%.6f\n", min_beta, max_beta);
    fprintf(fp, "input_parameters>\n");

    for (i = 0; i < N_REPLICAS; ++i)
    {
      betas[i] = min_beta + (max_beta - min_beta) * i / (N_REPLICAS - 1);
    }

    // 2. 初期化
    srand(48);
    int max_threads = omp_get_max_threads();
    pcg32_random_t rngs[max_threads];
    for (i = 0; i < max_threads; ++i)
    {
      pcg32_srandom_r(&rngs[i], time(NULL), (intptr_t)&rngs[i]);
    }
    initialize_replicas(replicas, N_REPLICAS, NX, NY, betas);

    // 全レプリカの物理量を格納する配列
    double total_Ex[N_REPLICAS] = {0.0};
    double total_Jx_squared[N_REPLICAS] = {0.0};
    double total_Ey[N_REPLICAS] = {0.0};
    double total_Jy_squared[N_REPLICAS] = {0.0};
    long measurement_count = 0;

    // 3. シミュレーションループ
    printf("<simulation_progress\n");
    sweep = 0;
    while (measurement_count < N_SAMPLING)
    {
      // 各レプリカでメトロポリス更新
  #pragma omp parallel for
      for (i = 0; i < N_REPLICAS; ++i)
      {
        metropolis_sweep(&replicas[i], &rngs[omp_get_thread_num()]);
      }

      // レプリカ交換
      for (i = 0; i < N_REPLICAS; ++i)
      { // 交換頻度を上げるためにループ
        replica_exchange(replicas, N_REPLICAS);
      }

      // 物理量の測定(初期緩和後)
      if (sweep > THERMALIZATION_SWEEPS && (sweep - THERMALIZATION_SWEEPS) % SAMPLE_INTERVAL_SWEEPS == 0)
      {
  #pragma omp parallel for
        for (i = 0; i < N_REPLICAS; ++i)
        {
          float E_x_sweep = 0.0f;
          float J_x_sweep = 0.0f;
          float E_y_sweep = 0.0f;
          float J_y_sweep = 0.0f;

  #pragma omp simd reduction(+ : E_x_sweep, J_x_sweep, E_y_sweep, J_y_sweep)
          for (int j = 0; j < N_SITES; ++j)
          {
            int ix = j % NX;
            int iy = j / NX;

            // x方向 (右の隣人)
            int right_neighbor_idx = iy * NX + (ix + 1) % NX;
            float delta_theta_x = replicas[i].spin[j] - replicas[i].spin[right_neighbor_idx];
            E_x_sweep += J * cosf(delta_theta_x);
            J_x_sweep += J * sinf(delta_theta_x);

            // y方向 (下の隣人)
            int down_neighbor_idx = ((iy + 1) % NY) * NX + ix;
            float delta_theta_y = replicas[i].spin[j] - replicas[i].spin[down_neighbor_idx];
            E_y_sweep += J * cosf(delta_theta_y);
            J_y_sweep += J * sinf(delta_theta_y);
          }
          total_Ex[i] += E_x_sweep;
          total_Jx_squared[i] += (double)J_x_sweep * J_x_sweep;
          total_Ey[i] += E_y_sweep;
          total_Jy_squared[i] += (double)J_y_sweep * J_y_sweep;
        }
        measurement_count++;
        if (measurement_count % (N_SAMPLING / 10) == 0)
        {
          printf("sweep=%d measurement_count=%ld total_Ex=%f total_Jx_squared=%f total_Ey=%f total_Jy_squared=%f\n",
                 sweep, measurement_count, total_Ex[0], total_Jx_squared[0], total_Ey[0], total_Jy_squared[0]);
        }
      }

      sweep++;
    }
    printf("simulation_progress>\n");

    // 最終結果の計算と表示
    fprintf(fp, "<helicity_modulus_results\n");
    for (i = 0; i < N_REPLICAS; ++i)
    {
      if (measurement_count > 0)
      {
        double avg_Ex = total_Ex[i] / measurement_count;
        double avg_Jx_squared = total_Jx_squared[i] / measurement_count;
        double upsilon_x = (avg_Ex / N_SITES) - (replicas[i].beta / N_SITES) * avg_Jx_squared;

        double avg_Ey = total_Ey[i] / measurement_count;
        double avg_Jy_squared = total_Jy_squared[i] / measurement_count;
        double upsilon_y = (avg_Ey / N_SITES) - (replicas[i].beta / N_SITES) * avg_Jy_squared;

        double helicity_modulus = (upsilon_x + upsilon_y) / 2.0;
        double energy_per_site = replicas[i].energy / (double)N_SITES;
        fprintf(fp, "replica=%d beta=%.6f energy_per_site=%.6f helicity_modulus=%.6f\n",
                i, replicas[i].beta, energy_per_site, helicity_modulus);
      }
    }
    fprintf(fp, "helicity_modulus_results>\n");

    // 4. メモリ解放
    free_replicas(replicas, N_REPLICAS);

    // end time measurement
    gettimeofday(&end_time, NULL);
    wall_time = (end_time.tv_sec - start_time.tv_sec) +
                (end_time.tv_usec - start_time.tv_usec) / 1000000.0;
    fprintf(fp, "<execution_time\n%.4f\nexecution_time>\n", wall_time);

    fclose(fp);

    return 0;
  }
  ```
]
