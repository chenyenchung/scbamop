package dedup

import "testing"

func TestMarkDuplicates(t *testing.T) {
	decisions := []readDecision{
		{
			ReadIndex:   0,
			CellBarcode: "CB1",
			UMI:         "UB1",
			ReferenceID: 0,
			Coordinate:  10,
			Strand:      0,
			MapQ:        20,
			Keep:        true,
		},
		{
			ReadIndex:   1,
			CellBarcode: "CB1",
			UMI:         "UB1",
			ReferenceID: 0,
			Coordinate:  10,
			Strand:      0,
			MapQ:        10,
			Keep:        true,
		},
		{
			ReadIndex:   2,
			CellBarcode: "CB1",
			UMI:         "UB2",
			ReferenceID: 0,
			Coordinate:  10,
			Strand:      0,
			MapQ:        30,
			Keep:        true,
		},
	}

	markDuplicates(decisions)

	if !decisions[0].Keep {
		t.Fatalf("expected first read to be kept")
	}
	if decisions[1].Keep {
		t.Fatalf("expected duplicate read to be discarded")
	}
	if !decisions[2].Keep {
		t.Fatalf("expected distinct molecule to be kept")
	}
}

func TestMarkDuplicatesTieBreak(t *testing.T) {
	decisions := []readDecision{
		{
			ReadIndex:   5,
			CellBarcode: "CB1",
			UMI:         "UB1",
			ReferenceID: 1,
			Coordinate:  200,
			Strand:      1,
			MapQ:        25,
			Keep:        true,
		},
		{
			ReadIndex:   2,
			CellBarcode: "CB1",
			UMI:         "UB1",
			ReferenceID: 1,
			Coordinate:  200,
			Strand:      1,
			MapQ:        25,
			Keep:        true,
		},
	}

	markDuplicates(decisions)

	if !decisions[0].Keep {
		t.Fatalf("expected earliest read index to be kept")
	}
	if decisions[1].Keep {
		t.Fatalf("expected later read index to be discarded")
	}
}
