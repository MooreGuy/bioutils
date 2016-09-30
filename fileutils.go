package bioutils

import (
	"bufio"
	"log"
	"os"
)

func ParseLines(filepath string) ([]string, error) {
	file, err := os.Open(filepath)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	var lines []string = make([]string, 5)
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		lines = append(lines, scanner.Text())
	}

	err = scanner.Err()
	return lines, err
}
