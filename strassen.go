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

func add(A, B [][]int) [][]int {
	n := len(A)
	m := len(A[0])

	C := make([][]int, n)
	for i := range C {
		C[i] = make([]int, m)
	}

	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			addCount++
			C[i][j] = A[i][j] + B[i][j]
		}
	}

	return C
}

func subtract(A, B [][]int) [][]int {
	n := len(A)
	m := len(A[0])

	C := make([][]int, n)
	for i := range C {
		C[i] = make([]int, m)
	}

	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			addCount++
			C[i][j] = A[i][j] - B[i][j]
		}
	}

	return C
}

func gauss(A [][]int, B [][]int) [][]int {
	n := len(A)
	m := len(B[0])

	C := make([][]int, n)
	for i := range C {
		C[i] = make([]int, m)
	}

	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			mulCount++
			C[i][j] = A[i][0] * B[0][j]
			for k := 1; k < len(A[0]); k++ {
				addCount++
				mulCount++
				C[i][j] += A[i][k] * B[k][j]
			}
		}
	}

	return C
}

func strassen(A, B [][]int) [][]int {
	n := len(A)

	if n == 1 {
		mulCount++
		return [][]int{{A[0][0] * B[0][0]}}
	}

	if n%2 == 1 {
		A11, A12, A21, A22 := splitIntoBlockMatricesDynamicPeeling(A)
		B11, B12, B21, B22 := splitIntoBlockMatricesDynamicPeeling(B)

		C11 := add(strassen(A11, B11), gauss(A12, B21))
		C12 := add(gauss(A11, B12), gauss(A12, B22))
		C21 := add(gauss(A21, B11), gauss(A22, B21))
		C22 := add(gauss(A21, B12), gauss(A22, B22))

		return concatenateMatrices(C11, C12, C21, C22)
	}

	A11, A12, A21, A22 := splitIntoBlockMatrices(A)
	B11, B12, B21, B22 := splitIntoBlockMatrices(B)

	P1 := strassen(add(A11, A22), add(B11, B22))
	P2 := strassen(add(A21, A22), B11)
	P3 := strassen(A11, subtract(B12, B22))
	P4 := strassen(A22, subtract(B21, B11))
	P5 := strassen(add(A11, A12), B22)
	P6 := strassen(subtract(A21, A11), add(B11, B12))
	P7 := strassen(subtract(A12, A22), add(B21, B22))

	C11 := add(subtract(add(P1, P4), P5), P7)
	C12 := add(P3, P5)
	C21 := add(P2, P4)
	C22 := add(subtract(add(P1, P3), P2), P6)

	return concatenateMatrices(C11, C12, C21, C22)
}

func runTimeTest(n int) []TestData {
	data := make([]TestData, 0)

	for i := 1; i < n; i++ {
		addCount = 0
		mulCount = 0
		A := generateRandomMatrix(i)
		B := generateRandomMatrix(i)

		t0 := time.Now()
		strassen(A, B)
		executionTime := time.Since(t0)
		fmt.Printf("Execution time for n = %v: %v\n", i, executionTime)
		flops := addCount + mulCount
		data = append(data, TestData{size: i, time: executionTime.Seconds(), flops: flops, add: addCount, mul: mulCount})
	}

	return data
}

func runAlgorithmTest(n int) {
	A := generateRandomMatrix(n)
	B := generateRandomMatrix(n)

	expected := gaussSkipCount(A, B)

	t0 := time.Now()
	actual := strassen(A, B)
	executionTime := time.Since(t0)
	fmt.Printf("Execution time for n = %v: %v\n", n, executionTime)

	assertMatrixMultiplicationIsCorrect(expected, actual)
	printInfo()
}

func main() {
	// runAlgorithmTest(511)
	data := runTimeTest(20)
	writeDataToFile(data)
}
