package main

import (
	"fmt"
	"time"
)

var mulCount, addCount int = 0, 0

type TestData struct {
	size  int
	time  float64
	flops int
	add   int
	mul   int
}

func add[T int | float64](A, B, C [][]T, xA1, xA2, yA1, yA2, xB1, xB2, yB1, yB2, xC1, xC2, yC1, yC2 int) {
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

func subtract[T int | float64](A, B, C [][]T, xA1, xA2, yA1, yA2, xB1, xB2, yB1, yB2, xC1, xC2, yC1, yC2 int) {
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

func gauss[T int | float64](A, B, C [][]T, xA1, xA2, yA1, yA2, xB1, xB2, yB1, yB2, xC1, xC2, yC1, yC2 int) {
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

func Strassen(A, B [][]float64) (C [][]float64) {
	n := len(A)
	C, tmp := initMatricesFloat64(n, 20)
	strassenRec(A, B, C, tmp, 0, n-1, 0, n-1, n, 0)
	return
}

func strassenRec[T int | float64](A, B, C [][]T, tmp [][][]T, x1, x2, y1, y2, size, depth int) {
	if size <= 5 {
		gauss(A, B, C, x1, x2, y1, y2, x1, x2, y1, y2, x1, x2, y1, y2)
		return
	}

	if size%2 == 1 {
		// Dynamic peeling split into block matrices
		// 11: x1, x2-1, y1, y2-1
		// 12: x1, x2-1, y2, y2
		// 21: x2, x2, y1, y2-1
		// 22: x2, x2, y2, y2

		// C11
		strassenRec(A, B, tmp[9*depth], tmp, x1, x2-1, y1, y2-1, size-1, depth+1)
		gauss(A, B, tmp[9*depth+1], x1, x2-1, y2, y2, x2, x2, y1, y2-1, x1, x2-1, y1, y2-1)
		add(tmp[9*depth], tmp[9*depth+1], C, x1, x2-1, y1, y2-1, x1, x2-1, y1, y2-1, x1, x2-1, y1, y2-1)

		// C12
		gauss(A, B, tmp[9*depth], x1, x2-1, y1, y2-1, x1, x2-1, y2, y2, x1, x2-1, y2, y2)
		gauss(A, B, tmp[9*depth+1], x1, x2-1, y2, y2, x2, x2, y2, y2, x1, x2-1, y2, y2)
		add(tmp[9*depth], tmp[9*depth+1], C, x1, x2-1, y2, y2, x1, x2-1, y2, y2, x1, x2-1, y2, y2)

		// C21
		gauss(A, B, tmp[9*depth], x2, x2, y1, y2-1, x1, x2-1, y1, y2-1, x2, x2, y1, y2-1)
		gauss(A, B, tmp[9*depth+1], x2, x2, y2, y2, x2, x2, y1, y2-1, x2, x2, y1, y2-1)
		add(tmp[9*depth], tmp[9*depth+1], C, x2, x2, y1, y2-1, x2, x2, y1, y2-1, x2, x2, y1, y2-1)

		// C22
		gauss(A, B, tmp[9*depth], x2, x2, y1, y2-1, x1, x2-1, y2, y2, x2, x2, y2, y2)
		gauss(A, B, tmp[9*depth+1], x2, x2, y2, y2, x2, x2, y2, y2, x2, x2, y2, y2)
		add(tmp[9*depth], tmp[9*depth+1], C, x2, x2, y2, y2, x2, x2, y2, y2, x2, x2, y2, y2)
		return
	}

	size_half := size / 2
	// Even split into block matrices
	// 11: x1, x1+size_half-1, y1, y1+size_half-1
	// 12: x1, x1+size_half-1, y1+size_half, y2
	// 21: x1+size_half, x2, y1, y1+size_half-1
	// 22: x1+size_half, x2, y1+size_half, y2

	// P1
	add(A, A, tmp[9*depth], x1, x1+size_half-1, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	add(B, B, tmp[9*depth+1], x1, x1+size_half-1, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	strassenRec(tmp[9*depth], tmp[9*depth+1], tmp[9*depth+2], tmp, x1, x1+size_half-1, y1, y1+size_half-1, size_half, depth+1)

	// P2
	add(A, A, tmp[9*depth], x1+size_half, x2, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	strassenRec(tmp[9*depth], B, tmp[9*depth+3], tmp, x1, x1+size_half-1, y1, y1+size_half-1, size_half, depth+1)

	// P3
	subtract(B, B, tmp[9*depth], x1, x1+size_half-1, y1+size_half, y2, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	strassenRec(A, tmp[9*depth], tmp[9*depth+4], tmp, x1, x1+size_half-1, y1, y1+size_half-1, size_half, depth+1)

	// P4
	subtract(B, B, tmp[9*depth], x1+size_half, x2, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2)
	strassenRec(A, tmp[9*depth], tmp[9*depth+5], tmp, x1+size_half, x2, y1+size_half, y2, size_half, depth+1)

	// P5
	add(A, A, tmp[9*depth], x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1+size_half, y2, x1+size_half, x2, y1+size_half, y2)
	strassenRec(tmp[9*depth], B, tmp[9*depth+6], tmp, x1+size_half, x2, y1+size_half, y2, size_half, depth+1)

	// P6
	subtract(A, A, tmp[9*depth], x1+size_half, x2, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1)
	add(B, B, tmp[9*depth+1], x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	strassenRec(tmp[9*depth], tmp[9*depth+1], tmp[9*depth+7], tmp, x1, x1+size_half-1, y1, y1+size_half-1, size_half, depth+1)

	// P7
	subtract(A, A, tmp[9*depth], x1, x1+size_half-1, y1+size_half, y2, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	add(B, B, tmp[9*depth+1], x1+size_half, x2, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	strassenRec(tmp[9*depth], tmp[9*depth+1], tmp[9*depth+8], tmp, x1, x1+size_half-1, y1, y1+size_half-1, size_half, depth+1)

	// C11
	add(tmp[9*depth+2], tmp[9*depth+5], tmp[9*depth], x1, x1+size_half-1, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	subtract(tmp[9*depth], tmp[9*depth+6], tmp[9*depth+1], x1, x1+size_half-1, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1, y1+size_half-1)
	add(tmp[9*depth+1], tmp[9*depth+8], C, x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1)

	// C12
	add(tmp[9*depth+4], tmp[9*depth+6], C, x1, x1+size_half-1, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2, x1, x1+size_half-1, y1+size_half, y2)

	// C21
	add(tmp[9*depth+3], tmp[9*depth+5], C, x1, x1+size_half-1, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2, x1+size_half, x2, y1, y1+size_half-1)

	// C22
	add(tmp[9*depth+2], tmp[9*depth+4], tmp[9*depth], x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1)
	subtract(tmp[9*depth], tmp[9*depth+3], tmp[9*depth+1], x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1)
	add(tmp[9*depth+1], tmp[9*depth+7], C, x1, x1+size_half-1, y1, y1+size_half-1, x1, x1+size_half-1, y1, y1+size_half-1, x1+size_half, x2, y1+size_half, y2)
}

func initMatricesInt(n, depth int) (result [][]int, tmp [][][]int) {
	result = make([][]int, n)
	for i := range result {
		result[i] = make([]int, n)
	}

	// hard-coded 9*depth: tmp1, tmp2, P1, ..., P7,
	tmp = make([][][]int, 9*depth)
	for i := 0; i < 9*depth; i++ {
		tmp[i] = make([][]int, n)
		for j := range tmp[i] {
			tmp[i][j] = make([]int, n)
		}
	}
	return
}

func initMatricesFloat64(n, depth int) (result [][]float64, tmp [][][]float64) {
	result = make([][]float64, n)
	for i := range result {
		result[i] = make([]float64, n)
	}

	// hard-coded 9*depth: tmp1, tmp2, P1, ..., P7
	tmp = make([][][]float64, 9*depth)
	for i := 0; i < 9*depth; i++ {
		tmp[i] = make([][]float64, n)
		for j := range tmp[i] {
			tmp[i][j] = make([]float64, n)
		}
	}
	return
}

func runTimeTest(n int) []TestData {
	data := make([]TestData, 0)

	for i := 1; i < n; i++ {
		addCount = 0
		mulCount = 0
		A := generateRandomMatrixFloat64(i)
		B := generateRandomMatrixFloat64(i)

		t0 := time.Now()
		Strassen(A, B)
		executionTime := time.Since(t0)
		fmt.Printf("Execution time for n = %v: %v\n", i, executionTime)
		flops := addCount + mulCount
		data = append(data, TestData{size: i, time: executionTime.Seconds(), flops: flops, add: addCount, mul: mulCount})
	}

	return data
}

func runAlgorithmTest(n int) {
	A := generateRandomMatrixFloat64(n)
	B := generateRandomMatrixFloat64(n)

	expected := gaussSkipCount(A, B)

	t0 := time.Now()
	actual := Strassen(A, B)
	executionTime := time.Since(t0)
	fmt.Printf("Execution time for n = %v: %v\n", n, executionTime)

	assertMatrixMultiplicationIsCorrect(expected, actual)
	printInfo()
}

func main() {
	// runAlgorithmTest(1500)
	data := runTimeTest(1501)
	writeDataToFile(data)
}
