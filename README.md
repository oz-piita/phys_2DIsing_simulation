# phys_2DIsing_simulation

二次元イジング模型のシミュレーションをGoで実装してみた。

実装はしたが正しくない。


# DEMO

工事中

# Features

マルコフニコフ連鎖モンテカルロ法に基づいたメトロポリス法と熱浴法で2次元イジングモデルをシミュレーションする。

さらに解析解を計算して、それぞれをプロットしてみる。

# Requirement

* Go 1.15.2
* GoogleColaboratry

# Usage


```bash
git clone https://github.com/oz-piita/phys_2DIsing_simulation/new/master
cd phys_2DIsing_simulation
```

シミュレーションを実行する
```bash
go run ./meropolis/mcmc_metro.go
go run ./gibbs/mcmc_gibbs.go
go run ./analitics/analitical_solution.go
```

Google colaboratryでgraph/plot_ising.ipynbを開き、Path変数を変更して実行

# Note

明らかに実装にミスが有ると思われる。直し方がわかりしだい直す。
* 比熱の標準偏差が大きすぎる
* 各シュミレーションで導かれる物理量のオーダーが100倍以上違っている


# Author

oz-piita
