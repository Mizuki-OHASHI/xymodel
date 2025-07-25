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

== 問題設定


