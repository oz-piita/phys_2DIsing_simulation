package main

// 2次元イジングモデルの解析解

import (
	"fmt"
	"math"
	"strconv"

	"encoding/csv"
	"log"
	"os"
)

const (
	j  = 1.0            // スピン結合定数
	kb = 1.3806 * 1e-23 // ボルツマン定数 1.3806 * 1e-23(J/K)
)

// 逆温度10点
var betas []float64 = []float64{
	0.11,
	0.2,
	0.64,
	1.17,
	1.7,
	2.23,
	2.76,
	3.29,
	3.82,
	4.88,
}

func main() {
	var ene, cap, mag []float64
	for _, beta := range betas {
		ene = append(ene, solveEnergy(beta))
		cap = append(cap, solveCapacity(beta))
		mag = append(mag, solveMag(beta))
	}
	tc := solveTC()
	fmt.Printf("転移温度Tc:%v\n", tc)
	fmt.Printf("内部エネルギー:%v\n", ene)
	fmt.Printf("比熱:%v\n", cap)
	fmt.Printf("自発磁化:%v\n", mag)
	// plot用にcsvを作成する
	outputCSV(sortTo2d(tc, ene, cap, mag))
}

// 第1種完全楕円積分
// complete elliptic Integral of the 1st kind
func ellip1st(kk float64) float64 {
	kap := (1 - math.Sqrt(1-kk*kk)) / (1 + math.Sqrt(1-kk*kk))
	coef := math.Pi * 0.5 * (1 + kap)

	temp := 1.0
	for n := 1; n < 10; n++ {
		mid := doubleFactrial(2*n-1) / doubleFactrial(2*n)
		temp += mid * mid * math.Pow(kap, float64(2*n))
	}
	return coef * temp
}

// 第2種完全楕円積分
// complete elliptic Integral of the 2nd kind
func ellip2nd(kk float64) float64 {
	kap := (1 - math.Sqrt(1-kk*kk)) / (1 + math.Sqrt(1-kk*kk))
	coef := math.Pi * 0.5 / (1 + kap)

	temp := 1.0
	for n := 1; n < 10; n++ {
		mid := doubleFactrial(2*n-3) / doubleFactrial(2*n)
		temp += mid * mid * math.Pow(kap, float64(2*n))
	}
	return coef * temp
}

// 二重階乗
func doubleFactrial(n int) float64 {
	m := 1
	for i := n; 0 < i; i -= 2 {
		m *= i
	}
	return float64(m)
}

// 内部エネルギー
func solveEnergy(beta float64) float64 {
	coef := 1 + 2*calcKddash(beta)*ellip1st(calcK(beta))/math.Pi
	return -j * coef / math.Tanh(2*j*beta)
}

// 比熱
func solveCapacity(beta float64) float64 {
	kk := calcK(beta)
	kkd := calcKddash(beta)
	coef1 := j * j * beta * beta * kb * 2 / math.Pi
	coef2invert := math.Tanh(2 * j * beta)
	form := 2*ellip1st(kk) - 2*ellip2nd(kk) - (1-kkd)*(math.Pi/2+kkd*ellip1st(kk))
	return coef1 * form / (coef2invert * coef2invert)
}

// 自発磁化
func solveMag(beta float64) float64 {
	tc := solveTC()
	if beta <= 1/(kb*tc) {
		return 0.0
	}
	sinh4 := math.Pow(math.Sinh(2*j*beta), -4)
	return math.Pow(1-sinh4, 1/8)
}

// 転移温度
func solveTC() float64 {
	return 2 * j / kb / math.Log(1+math.Sqrt(2))
}

// return k
func calcK(beta float64) float64 {
	d := 2 * j * beta
	return 2 * math.Tanh(d) / math.Cosh(d)
}

// return k"
func calcKddash(beta float64) float64 {
	d := math.Tanh(2 * j * beta)
	return 2*d*d - 1
}

func sortTo2d(tc float64, ene, cap, mag []float64) [][]string {
	records := make([][]string, 11)
	records[0] = []string{"逆温度", "内部エネルギー", "比熱", "自発磁化", "転移温度"}
	for i, _ := range betas {
		r := make([]string, 5)
		r[0] = strconv.FormatFloat(betas[i], 'f', -1, 64)
		r[1] = strconv.FormatFloat(ene[i], 'f', -1, 64)
		r[2] = strconv.FormatFloat(cap[i], 'f', -1, 64)
		r[3] = strconv.FormatFloat(mag[i], 'f', -1, 64)
		r[4] = strconv.FormatFloat(tc, 'f', -1, 64)
		records[i+1] = r
	}
	return records
}

func outputCSV(records [][]string) {

	f, err := os.Create("./result/analitic.csv")
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
