package metadata

import "testing"

func TestSanitizeLabel(t *testing.T) {
	sanitized, modified := sanitizeLabel("bad/label")
	if sanitized != "bad_label" {
		t.Fatalf("expected sanitized label, got %q", sanitized)
	}
	if !modified {
		t.Fatalf("expected modification flag")
	}

	sanitized, modified = sanitizeLabel(".hidden")
	if sanitized != "_hidden" {
		t.Fatalf("expected leading dot replacement, got %q", sanitized)
	}
	if !modified {
		t.Fatalf("expected modification flag")
	}

	sanitized, modified = sanitizeLabel("foo..bar")
	if sanitized != "foo__bar" {
		t.Fatalf("expected '..' replacement, got %q", sanitized)
	}
	if !modified {
		t.Fatalf("expected modification flag")
	}
}
