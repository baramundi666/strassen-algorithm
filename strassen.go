package main

import (
	"fmt"
	"time"
)

var mulCount, addCount int = 0, 0

var result [][]int
var P1, P2, P3, P4, P5, P6, P7 [][][]int
var tmp1, tmp2, tmp3 [][][]int

type TestData struct {
	size  int
	time  float64
	flops int
	add   int
	mul   int
}

func add(A, B, C [][]int, xA1, xA2, yA1, yA2, xB1, xB2, yB1, yB2, xC1, xC2, yC1, yC2 int) {
	// fmt.Println("add")
	// printMatrix(A, xA1, xA2, yA1, yA2)
	// printMatrix(B, xB1, xB2, yB1, yB2)
	// fmt.Println(xC1, xC2, yC1, yC2)
	// fmt.Println()
	for i := 0; i <= xA2-xA1; i++ {
		for j := 0; j <= yA2-yA1; j++ {
			addCount++
			C[xC1+i][yC1+j] = A[xA1+i][yA1+j] + B[xB1+i][yB1+j]
		}
	}
}

func subtract(A, B, C [][]int, xA1, xA2, yA1, yA2, xB1, xB2, yB1, yB2, xC1, xC2, yC1, yC2 int) {
	// fmt.Println("sub")
	// printMatrix(A, xA1, xA2, yA1, yA2)
	// printMatrix(B, xB1, xB2, yB1, yB2)
	// fmt.Println(xC1, xC2, yC1, yC2)
	// fmt.Println()
	for i := 0; i <= xA2-xA1; i++ {
		for j := 0; j <= yA2-yA1; j++ {
			addCount++
			C[xC1+i][yC1+j] = A[xA1+i][yA1+j] - B[xB1+i][yB1+j]
		}
	}
}

func gauss(A, B, C [][]int, xA1, xA2, yA1, yA2, xB1, xB2, yB1, yB2, xC1, xC2, yC1, yC2 int) {
	// fmt.Println("gauss")
	// printMatrix(A, xA1, xA2, yA1, yA2)
	// printMatrix(B, xB1, xB2, yB1, yB2)
	for i := 0; i <= xA2-xA1; i++ {
		for j := 0; j <= yB2-yB1; j++ {
			mulCount++
			C[xC1+i][yC1+j] = A[xA1+i][yA1] * B[yA1][yB1+j]
			for k := 1; k <= yA2-yA1; k++ {
				addCount++
				mulCount++
				C[xC1+i][yC1+j] += A[xA1+i][yA1+k] * B[yA1+k][yB1+j]
			}
		}
	}
}

func strassen(A, B, C [][]int, x1, x2, y1, y2, size, depth int) {
	if size <= 5 {
		gauss(A, B, C, x1, x2, y1, y2, x1, x2, y1, y2, x1, x2, y1, y2)
		return
	}
	// if size == 1 {
	// 	C[x1][y1] = A[x1][y1] * B[x1][y1]
	// 	return
	// }

	if size%2 == 1 {
		// Dynamic peeling split into block matrices
		// 11: x1, x2-1, y1, y2-1
		// 12: x1, x2-1, y2, y2
		// 21: x2, x2, y1, y2-1
		// 22: x2, x2, y2, y2

		// C11
		strassen(A, B, tmp3[depth], x1, x2-1, y1, y2-1, size-1, depth+1)
		gauss(A, B, tmp2[depth], x1, x2-1, y2, y2, x2, x2, y1, y2-1, x1, x2-1, y1, y2-1)
		add(tmp3[depth], tmp2[depth], C, x1, x2-1, y1, y2-1, x1, x2-1, y1, y2-1, x1, x2-1, y1, y2-1)

		// C12
		gauss(A, B, tmp1[depth], x1, x2-1, y1, y2-1, x1, x2-1, y2, y2, x1, x2-1, y2, y2)
		gauss(A, B, tmp2[depth], x1, x2-1, y2, y2, x2, x2, y2, y2, x1, x2-1, y2, y2)
		add(tmp1[depth], tmp2[depth], C, x1, x2-1, y2, y2, x1, x2-1, y2, y2, x1, x2-1, y2, y2)

		// C21
		gauss(A, B, tmp1[depth], x2, x2, y1, y2-1, x1, x2-1, y1, y2-1, x2, x2, y1, y2-1)
		gauss(A, B, tmp2[depth], x2, x2, y2, y2, x2, x2, y1, y2-1, x2, x2, y1, y2-1)
		add(tmp1[depth], tmp2[depth], C, x2, x2, y1, y2-1, x2, x2, y1, y2-1, x2, x2, y1, y2-1)

		// C22
		gauss(A, B, tmp1[depth], x2, x2, y1, y2-1, x1, x2-1, y2, y2, x2, x2, y2, y2)
		gauss(A, B, tmp2[depth], x2, x2, y2, y2, x2, x2, y2, y2, x2, x2, y2, y2)
		add(tmp1[depth], tmp2[depth], C, x2, x2, y2, y2, x2, x2, y2, y2, x2, x2, y2, y2)
		return
	}

	size_half := size / 2
	// Even split into block matrices
	// 11: x1, x1+size_half-1, y1, y1+size_half-1
	// 12: x1, x1+size_half-1, y1+size_half, y2
	// 21: x1+size_half, x2, y1, y1+size_half-1
	// 22: x1+size_half, x2, y1+size_half, y2

	// P1
	add(A, A, tmp1[depth], x1, x1+size_half-1, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	add(B, B, tmp2[depth], x1, x1+size_half-1, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	strassen(tmp1[depth], tmp2[depth], P1[depth], x1, x1+size_half-1, y1, y1+size_half-1, size_half, depth+1)

	// P2
	add(A, A, tmp1[depth], x1+size_half, x2, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	strassen(tmp1[depth], B, P2[depth], x1, x1+size_half-1, y1, y1+size_half-1, size_half, depth+1)

	// P3
	subtract(B, B, tmp1[depth], x1, x1+size_half-1, y1+size_half, y2, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	strassen(A, tmp1[depth], P3[depth], x1, x1+size_half-1, y1, y1+size_half-1, size_half, depth+1)

	// P4
	subtract(B, B, tmp1[depth], x1+size_half, x2, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2)
	strassen(A, tmp1[depth], P4[depth], x1+size_half, x2, y1+size_half, y2, size_half, depth+1)

	// P5
	add(A, A, tmp1[depth], x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1+size_half, y2, x1+size_half, x2, y1+size_half, y2)
	strassen(tmp1[depth], B, P5[depth], x1+size_half, x2, y1+size_half, y2, size_half, depth+1)

	// P6
	subtract(A, A, tmp1[depth], x1+size_half, x2, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1)
	add(B, B, tmp2[depth], x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	strassen(tmp1[depth], tmp2[depth], P6[depth], x1, x1+size_half-1, y1, y1+size_half-1, size_half, depth+1)

	// P7
	subtract(A, A, tmp1[depth], x1, x1+size_half-1, y1+size_half, y2, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	add(B, B, tmp2[depth], x1+size_half, x2, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	strassen(tmp1[depth], tmp2[depth], P7[depth], x1, x1+size_half-1, y1, y1+size_half-1, size_half, depth+1)

	// C11
	add(P1[depth], P4[depth], tmp1[depth], x1, x1+size_half-1, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	subtract(tmp1[depth], P5[depth], tmp2[depth], x1, x1+size_half-1, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	add(tmp2[depth], P7[depth], C, x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1)

	// C12
	add(P3[depth], P5[depth], C, x1, x1+size_half-1, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1+size_half, y2)

	// C21
	add(P2[depth], P4[depth], C, x1, x1+size_half-1, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2, x1+size_half, x2, y1, y1+size_half-1)

	// C22
	add(P1[depth], P3[depth], tmp1[depth], x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1)
	subtract(tmp1[depth], P2[depth], tmp2[depth], x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1)
	add(tmp2[depth], P6[depth], C, x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2)
}

func initMatrices(n int) {
	result = make([][]int, n)
	for i := range result {
		result[i] = make([]int, n)
	}
	depth := 20
	P1 = make([][][]int, depth)
	P2 = make([][][]int, depth)
	P3 = make([][][]int, depth)
	P4 = make([][][]int, depth)
	P5 = make([][][]int, depth)
	P6 = make([][][]int, depth)
	P7 = make([][][]int, depth)
	tmp1 = make([][][]int, depth)
	tmp2 = make([][][]int, depth)
	tmp3 = make([][][]int, depth)

	for j := 0; j < 20; j++ {
		P1[j] = make([][]int, n)
		for i := range P1[j] {
			P1[j][i] = make([]int, n)
		}

		P2[j] = make([][]int, n)
		for i := range P2[j] {
			P2[j][i] = make([]int, n)
		}

		P3[j] = make([][]int, n)
		for i := range P3[j] {
			P3[j][i] = make([]int, n)
		}

		P4[j] = make([][]int, n)
		for i := range P4[j] {
			P4[j][i] = make([]int, n)
		}

		P5[j] = make([][]int, n)
		for i := range P5[j] {
			P5[j][i] = make([]int, n)
		}

		P6[j] = make([][]int, n)
		for i := range P6[j] {
			P6[j][i] = make([]int, n)
		}

		P7[j] = make([][]int, n)
		for i := range P7[j] {
			P7[j][i] = make([]int, n)
		}

		tmp1[j] = make([][]int, n)
		for i := range tmp1[j] {
			tmp1[j][i] = make([]int, n)
		}

		tmp2[j] = make([][]int, n)
		for i := range tmp2[j] {
			tmp2[j][i] = make([]int, n)
		}

		tmp3[j] = make([][]int, n)
		for i := range tmp3[j] {
			tmp3[j][i] = make([]int, n)
		}
	}
}

func runTimeTest(n int) []TestData {
	data := make([]TestData, 0)

	initMatrices(n)

	for i := 1; i < n; i++ {
		addCount = 0
		mulCount = 0
		A := generateRandomMatrix(i)
		B := generateRandomMatrix(i)

		t0 := time.Now()
		strassen(A, B, result, 0, i-1, 0, i-1, i, 0)
		executionTime := time.Since(t0)
		fmt.Printf("Execution time for n = %v: %v\n", i, executionTime)
		flops := addCount + mulCount
		data = append(data, TestData{size: i, time: executionTime.Seconds(), flops: flops, add: addCount, mul: mulCount})
	}

	return data
}

func runAlgorithmTest(n int) {
	initMatrices(n)

	A := generateRandomMatrix(n)
	B := generateRandomMatrix(n)

	expected := gaussSkipCount(A, B)

	t0 := time.Now()
	strassen(A, B, result, 0, n-1, 0, n-1, n, 0)
	executionTime := time.Since(t0)
	fmt.Printf("Execution time for n = %v: %v\n", n, executionTime)

	assertMatrixMultiplicationIsCorrect(expected, result)
	printInfo()
}

func main() {
	// runAlgorithmTest(1500)
	data := runTimeTest(1501)
	writeDataToFile(data)
}
