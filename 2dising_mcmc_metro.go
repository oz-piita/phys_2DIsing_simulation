package main

import (
	"fmt"
	"math"
	"math/rand"
	"time"
)

// グリッドNと外部磁場h
const Nx, Ny, h = 100, 100, 0.01

// 逆温度beta.(0.01, 5.0)の範囲で十点
var betas []float64 = []float64{0.2, 0.7, 1.2, 1.7, 2.2, 2.7, 3.2, 3.7, 4.2, 4.9}
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
	bMag, bMag_err, bEne, bEne_err, bCap := mcmcBeta(s)
	fmt.Println(bMag)
	fmt.Println(bMag_err)
	fmt.Println(bEne)
	fmt.Println(bEne_err)
	fmt.Println(bCap)
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
	return s[x_right][y] + s[x_left][y] + s[x][y_above] + s[x][y_bottom]
	// 近接スピン和
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

// metropolis法で掻き乱す
func metropolis(s [][]int, beta float64) [][]int {
	xShuffled, yShuffled := createShuffledInt(Nx), createShuffledInt(Ny)
	for x := range xShuffled {
		for y := range yShuffled {
			k := float64(neighborSpinSum(s, x, y)) + h
			trans_prob := math.Pow(math.E, float64(-2.0*beta*float64(s[x][y])*k))
			if rand.Float64() <= trans_prob {
				s[x][y] = -s[x][y]
			}
		}
	}
	return s
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
	interval := 10
	burn_in := 100
	for i := 0; i < 1000; i++ {
		s = metropolis(s, beta) // メトロポリス法で配位を確定
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

// 逆温度ごとににmcmcで物理量を取り出す
func mcmcBeta(s [][]int) (bMag, bMag_err, bEne, bEne_err, bCap []float64) {
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
		// bCap_err = append(bCap_err,errFloat())
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

// 標準”誤”差.エラーバー.
func errFloat(a []float64) float64 {
	return math.Sqrt(distFloat(a)) / math.Sqrt(float64(len(a)-1))
}
