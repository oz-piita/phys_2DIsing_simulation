package main

// 二次元イジング模型の熱浴法（ギブスサンプリング法）

import (
	"fmt"
	"math"
	"math/rand"
	"strconv"
	"time"

	"encoding/csv"
	"log"
	"os"
)

// グリッドNと外部磁場h
const Nx, Ny, h = 100, 100, 0.01

// 逆温度beta.(0.01, 5.0)の範囲で10点初期化
var betas []float64 = []float64{0.11, 0.2, 0.64, 1.17, 1.7, 2.23, 2.76, 3.29, 3.82, 4.88}
var s = make([][]int, Nx, Ny) // spin配位

func main() {
	// 初期配置を乱数で生成
	rand.Seed(time.Now().UnixNano())
	for i := 0; i < Nx; i++ {
		l := make([]int, Ny)
		for j := 0; j < Ny; j++ {
			l[j] = rand.Intn(2)
		}
		s[i] = l
	}
	// 順に磁化、磁化の標準偏差、エネルギー、その標準偏差、比熱
	bMag, bMag_err, bEne, bEne_err, bCap, bCap_err := mcmcBeta(s)
	fmt.Println(bMag, bMag_err)
	fmt.Println(bEne, bEne_err)
	fmt.Println(bCap, bCap_err)
	// plot用にcsvを作成する
	outputCSV(sortTo2d(bEne, bEne_err, bCap, bCap_err, bMag, bMag_err))
}

func neighborSpinSum(s [][]int, x int, y int) int {
	// 近接格子座標と周期境界条件
	x_right := x + 1
	if Nx <= x_right {
		x_right -= Nx
	}
	x_left := x - 1
	if x_left < 0 {
		x_left += Nx
	}
	y_above := y + 1
	if Ny <= y_above {
		y_above -= Ny
	}
	y_bottom := y - 1
	if y_bottom < 0 {
		y_bottom += Ny
	}
	// 近接スピン和
	return s[x_right][y] + s[x_left][y] + s[x][y_above] + s[x][y_bottom]
}

func calcEnergy(s [][]int) (energy float64) {
	for i := 0; i < Nx; i++ {
		for j := 0; j < Ny; j++ {
			energy += float64(-neighborSpinSum(s, i, j)) / 2
		}
	}
	sum := 0.0 // スピンの総和
	for _, v := range s {
		for _, w := range v {
			sum += float64(w)
		}
	}
	energy += h * sum // 外部磁場の効果
	return
}

// 熱浴法.gibbs sampling method.
func gibbs(s [][]int, beta float64) [][]int {
	xShuffled, yShuffled := createShuffledInt(Nx), createShuffledInt(Ny)
	for x := range xShuffled {
		for y := range yShuffled {
			k := float64(neighborSpinSum(s, x, y)) - h
			trans_prob := exp(beta*k) / (exp(beta*k) + exp(-beta*k))
			if rand.Float64() <= trans_prob {
				s[x][y] = 1
			} else {
				s[x][y] = -1
			}
		}
	}
	return s
}

func exp(a float64) float64 {
	return math.Pow(math.E, a)
}

// n個の連続整数の順序をランダムに入れ替えた1次元スライスを生成する
func createShuffledInt(n int) (x []int) {
	for i := 0; i < 100; i++ {
		x = append(x, i)
	}
	rand.Seed(time.Now().UnixNano())
	rand.Shuffle(len(x), func(i, j int) { x[i], x[j] = x[j], x[i] })
	return
}

// マルコフ連鎖モンテカルロ法
func mcmc(s [][]int, beta float64) (ms, energies, squareEne []float64) {
	steps := 1000  // モンテカルロstep
	interval := 10 // 実際にサンプルを採取するstep周期
	burn_in := 100 // サンプル採取前のギブス法の実行step数.バーンイン時間
	for i := 0; i < steps; i++ {
		s = gibbs(s, beta) // ギブス法で配位を確定
		if i%interval == 0 && burn_in <= i {
			m := 0.0 // 磁化
			for i := 0; i < Nx; i++ {
				for j := 0; j < Ny; j++ {
					m += float64(s[i][j])
				}
			}
			ms = append(ms, m)
			energy := calcEnergy(s)
			energies = append(energies, energy)
			squareEne = append(squareEne, energy*energy)
		}
	}
	return
}

// 逆温度ごとににmcmcで物理量とその標準誤差を取り出す
func mcmcBeta(s [][]int) (bMag, bMag_err, bEne, bEne_err, bCap, bCap_err []float64) {
	for _, beta := range betas {
		ms, energies, square_energies := mcmc(s, beta)

		bm := meanFloat(ms)
		bMag = append(bMag, bm)
		bMag_err = append(bMag_err, errFloat(ms))

		benergy := meanFloat(energies)
		bEne = append(bEne, benergy)
		bEne_err = append(bEne_err, errFloat(energies))

		bsqene := meanFloat(square_energies)
		bcapacity := (beta * beta) * (bsqene - benergy*benergy)
		bCap = append(bCap, bcapacity)

		bCap_err = append(bCap_err, beta*beta*errFloat(square_energies))
		// エネルギー平方からの誤差の伝播ってこれでいいのでしたっけ…自信がありません
	}
	return
}

// 和
func sumFloat(a []float64) (sum float64) {
	for _, v := range a {
		sum += v
	}
	return
}

// 平均
func meanFloat(a []float64) float64 {
	return sumFloat(a) / float64(len(a))
}

// 分散
func distFloat(a []float64) float64 {
	d := 0.0
	xm := meanFloat(a)
	for i := range a {
		d += (a[i] - xm) * (a[i] - xm)
	}
	return d / float64(len(a))
}

// 標準誤差.エラーバー.
func errFloat(a []float64) float64 {
	return math.Sqrt(distFloat(a)) / math.Sqrt(float64(len(a)-1))
}

func sortTo2d(ene, eneErr, cap, capErr, mag, magErr []float64) [][]string {
	records := make([][]string, 11)
	records[0] = []string{"逆温度", "内部エネルギー", "内部エネルギー標準誤差", "比熱", "比熱標準誤差", "自発磁化", "磁化標準誤差"}
	for i, _ := range betas {
		r := make([]string, 7)
		r[0] = strconv.FormatFloat(betas[i], 'f', -1, 64)
		r[1] = strconv.FormatFloat(ene[i], 'f', -1, 64)
		r[2] = strconv.FormatFloat(eneErr[i], 'f', -1, 64)
		r[3] = strconv.FormatFloat(cap[i], 'f', -1, 64)
		r[4] = strconv.FormatFloat(capErr[i], 'f', -1, 64)
		r[5] = strconv.FormatFloat(mag[i], 'f', -1, 64)
		r[6] = strconv.FormatFloat(magErr[i], 'f', -1, 64)
		records[i+1] = r
	}
	return records
}

func outputCSV(records [][]string) {
	f, err := os.Create("./result/gibbs.csv")
	if err != nil {
		log.Fatal(err)
	}
	w := csv.NewWriter(f)
	w.WriteAll(records)
	w.Flush()
	if err := w.Error(); err != nil {
		log.Fatal(err)
	}
}
