package main

import (
	"os"
	"os/exec"
	"testing"
)

func TestAll(t *testing.T) {
	if os.Getenv("SCBAMOP_TEST_ALL") != "" {
		return
	}

	cmd := exec.Command("go", "test", "./...")
	cmd.Env = append(os.Environ(), "SCBAMOP_TEST_ALL=1")
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("go test ./... failed: %v", err)
	}
}
