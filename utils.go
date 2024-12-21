package main

import (
	"fmt"
	"log"
	"math/rand"
	"os"
	"strconv"
	"time"
)

func splitIntoBlockMatrices(A [][]int) ([][]int, [][]int, [][]int, [][]int) {
	n_half := len(A) / 2

	A11, A12, A21, A22 := make([][]int, n_half), make([][]int, n_half), make([][]int, n_half), make([][]int, n_half)

	for i := 0; i < n_half; i++ {
		A11[i] = A[i][:n_half]
		A12[i] = A[i][n_half:]
		A21[i] = A[i+n_half][:n_half]
		A22[i] = A[i+n_half][n_half:]
	}

	return A11, A12, A21, A22
}

func splitIntoBlockMatricesDynamicPeeling(A [][]int) ([][]int, [][]int, [][]int, [][]int) {
	n := len(A)

	A11 := make([][]int, n-1)
	A12 := make([][]int, n-1)
	A21 := make([][]int, 1)
	A22 := [][]int{{A[n-1][n-1]}}

	for i := 0; i < n-1; i++ {
		A11[i] = A[i][:n-1]
		A12[i] = []int{A[i][n-1]}
	}

	A21[0] = A[n-1][:n-1]

	return A11, A12, A21, A22
}

func concatenateMatrices(A, B, C, D [][]int) [][]int {
	result := make([][]int, len(A)+len(C))

	for i := 0; i < len(A); i++ {
		result[i] = append(A[i], B[i]...)
	}

	for i := 0; i < len(C); i++ {
		result[len(A)+i] = append(C[i], D[i]...)
	}

	return result
}

func gaussSkipCount(A [][]int, B [][]int) [][]int {
	n := len(A)
	m := len(B[0])

	C := make([][]int, n)
	for i := range C {
		C[i] = make([]int, m)
	}

	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			C[i][j] = A[i][0] * B[0][j]
			for k := 1; k < len(A[0]); k++ {
				C[i][j] += A[i][k] * B[k][j]
			}
		}
	}

	return C
}

func generateRandomMatrix(n int) [][]int {
	A := make([][]int, n)
	for i := range A {
		A[i] = make([]int, n)
	}

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			A[i][j] = 1 + rand.Intn(7)
		}
	}

	return A
}

func assertMatrixMultiplicationIsCorrect(expected, actual [][]int) {
	n := len(expected)

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if expected[i][j] != actual[i][j] {
				log.Fatalf("Matrix multiplication error!\nexpected=%v\nactual=%v\n", expected[i][j], actual[i][j])
			}
		}
	}
	fmt.Println("Matrix multiplication was successful!")
}

func printMatrix(A [][]int) {
	for _, row := range A {
		fmt.Println(row)
	}
	fmt.Println()
}

func printInfo() {
	fmt.Printf("Multiplications: %v\n", mulCount)
	fmt.Printf("Additions: %v\n", addCount)
	flops := addCount + mulCount
	fmt.Printf("FLOPs: %v\n", flops)
}

func writeDataToFile(data []TestData) error {
	currentTime := time.Now()
	_, month, day := currentTime.Date()
	hour, minute, _ := currentTime.Clock()
	filename := "data/strassen-" + strconv.Itoa(day) + month.String() + "-" + strconv.Itoa(hour) + strconv.Itoa(minute) + ".csv"
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	if _, err := fmt.Fprintf(file, "size,time,flops,add,mul\n"); err != nil {
		return err
	}

	for _, testData := range data {
		if _, err := fmt.Fprintf(file, "%d,%f,%d,%d,%d\n", testData.size, testData.time, testData.flops, testData.add, testData.mul); err != nil {
			return err
		}
	}

	return nil
}
