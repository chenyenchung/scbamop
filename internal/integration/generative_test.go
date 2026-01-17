package integration

import (
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

func TestGenerativeSuite(t *testing.T) {
	if os.Getenv("RUN_GENERATIVE_TESTS") == "" {
		t.Skip("RUN_GENERATIVE_TESTS not set")
	}

	repoRoot, err := findRepoRoot()
	if err != nil {
		t.Fatalf("failed to locate repo root: %v", err)
	}

	pythonPath := os.Getenv("SCBAMOP_PYTHON")
	if pythonPath == "" {
		pythonPath = filepath.Join(repoRoot, ".venv", "bin", "python")
	}
	if _, err := os.Stat(pythonPath); err != nil {
		t.Skip("python interpreter not found for generative tests")
	}

	if err := checkPysam(pythonPath, repoRoot); err != nil {
		t.Skipf("pysam not available: %v", err)
	}

	binaryPath := filepath.Join(t.TempDir(), "scbamop")
	buildCmd := exec.Command("go", "build", "-o", binaryPath, ".")
	buildCmd.Dir = repoRoot
	if output, err := buildCmd.CombinedOutput(); err != nil {
		t.Fatalf("failed to build scbamop: %v\n%s", err, output)
	}

	scriptPath := filepath.Join(repoRoot, "tests", "scripts", "run_test.sh")
	modes := []string{"split", "dedup"}
	for _, mode := range modes {
		workDir := filepath.Join(t.TempDir(), mode)
		cmd := exec.Command("bash", scriptPath,
			"--python", pythonPath,
			"--scbamop", binaryPath,
			"--work-dir", workDir,
			"--mode", mode,
			"--unsafe-label",
		)
		cmd.Dir = repoRoot
		output, err := cmd.CombinedOutput()
		if err != nil {
			t.Fatalf("generative %s test failed: %v\n%s", mode, err, output)
		}
	}
}

func checkPysam(pythonPath string, repoRoot string) error {
	cmd := exec.Command(pythonPath, "-c", "import pysam")
	cmd.Dir = repoRoot
	return cmd.Run()
}

func findRepoRoot() (string, error) {
	current, err := os.Getwd()
	if err != nil {
		return "", err
	}
	for {
		if _, err := os.Stat(filepath.Join(current, "go.mod")); err == nil {
			return current, nil
		}
		parent := filepath.Dir(current)
		if parent == current {
			return "", os.ErrNotExist
		}
		current = parent
	}
}
