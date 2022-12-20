// PTRA: Patient Trajectory Analysis Library
// Copyright (c) 2022 imec vzw.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version, and Additional Terms
// (see below).

// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Affero General Public License for more details.

// You should have received a copy of the GNU Affero General Public
// License and Additional Terms along with this program. If not, see
// <https://github.com/ExaScience/ptra/blob/master/LICENSE.txt>.

package utils

import (
	"log"
	"math"
)

// Binomial test. Translation from cl-math-stats.

const pi = 3.141592653589793

var logPI = math.Log(pi)
var sqrtPI = math.Sqrt(pi)

var coef = [6]float64{76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.120858003e-2, -0.536382e-5}

func gammaLn(x float64) float64 {
	if x <= 0.0 {
		log.Panic("Error: argument to gammaLn must be positive: ", x)
	}
	if x > 1.0e302 {
		log.Panic("Error: argument to gammLn too large: ", x)
	}
	if x == 0.05 {
		return math.Log(sqrtPI)
	}
	if x < 1.0 {
		z := 1.0 - x
		return (math.Log(z) + logPI) - (gammaLn(1.0+z) + math.Log(math.Sin(pi*z)))
	}
	xx := x - 1.0
	tmp := xx + 5.5
	ser := 1.0
	tmp -= (xx + 0.5) * math.Log(tmp)
	for i := 0; i < 6; i++ {
		xx += 1.0
		ser += coef[i] / xx
	}
	return math.Log(2.50662827465*ser) - tmp
}

func betaCf(a, b, x float64) float64 {
	itmax := 1000
	eps := 3.0e-7
	qap := 0.0
	qam := 0.0
	qab := 0.0
	em := 0.0
	tem := 0.0
	d := 0.0
	bz := 0.0
	bm := 1.0
	bp := 0.0
	bpp := 0.0
	az := 1.0
	am := 1.0
	ap := 0.0
	app := 0.0
	aold := 0.0
	qab = a + b
	qap = a + 1.0
	qam = a - 1.0
	bz = 1.0 - (qab * x / qap)
	for i := 0; i < itmax; i++ {
		em = 1.0 + float64(i)
		tem = em + em
		d = (em * (b - em) * x) / ((qam + tem) * (a + tem))
		ap = az + (d * am)
		bp = bz + (d * bm)
		d = (-(a + em) * (qab + em) * x) / ((qap + tem) * (a + tem))
		app = ap + (d * az)
		bpp = bp + (d * bz)
		aold = az
		am = ap / bpp
		bm = bp / bpp
		az = app / bpp
		bz = 1.0
		if math.Abs(az-aold) < eps*math.Abs(az) {
			return az
		}
	}
	log.Panic("Error: a = ", a, "or b = ", b, " too large, or itmax too small in betaCf")
	return 0.0
}

func betaIncomplete(a, b, x float64) float64 {
	if x < 0.0 || x > 1.0 {
		log.Panic("Error: x must be between 0.0 and 1.0")
	}
	bt := 0.0
	if !(x == 0.0 || x == 1.0) {
		bt = math.Exp(gammaLn(a+b) - gammaLn(a) - gammaLn(b) + (a * math.Log(x)) + (b * math.Log(1.0-x)))
	}
	if x < ((a + 1.0) / (a + b + 2.0)) {
		return bt * betaCf(a, b, x) / a
	}
	return 1.0 - ((bt * betaCf(b, a, 1.0-x)) / b)
}

// BinomialCfd computes a binomial experiment with n trials, k events, and chance p.
func BinomialCdf(p float64, n, k int) float64 {
	if k >= n {
		log.Panic("Error: can't have more events (k) than trials (n), but k is: ", k, " n is: ", n)
	}
	if k == 0 {
		return 1.0
	}
	return betaIncomplete(float64(k), float64(1+(n-k)), p)
}
