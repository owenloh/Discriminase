"""2-bit DNA encoding.

A genome is stored as a ``bitarray`` with two bits per base (A=00, C=01, G=10,
T=11). The reverse complement is precomputed once so both strands can be sliced
cheaply. This packing is what makes the Hamming distance an XOR + popcount and
cuts memory ~4x versus an ASCII string.
"""
from bitarray import bitarray

from .config import GuideFinderConfig


class TwoBitDNA:
    config = GuideFinderConfig()
    if config.all_N_into_A:
        # Ambiguous N bases are folded into A. This is a deliberate simplification:
        # it lets every base encode in 2 bits, at the cost of treating N as A.
        _encode_map = {'A': bitarray('00'), 'C': bitarray('01'), 'G': bitarray('10'), 'T': bitarray('11'), 'N': bitarray('00')}
        _decode_map = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
        _comp_encode_map = {'A': bitarray('11'), 'C': bitarray('10'), 'G': bitarray('01'), 'T': bitarray('00'), 'N': bitarray('11')}
    else:
        _encode_map = {'A': bitarray('00'), 'C': bitarray('01'), 'G': bitarray('10'), 'T': bitarray('11')}
        _decode_map = {v.to01(): k for k, v in _encode_map.items()}
        _comp_encode_map = {'A': bitarray('11'), 'C': bitarray('10'), 'G': bitarray('01'), 'T': bitarray('00')}

    organism_count = 0

    def __init__(self, name: str, sequence: str):
        self.name = name
        self.organism_num = self.__class__.organism_count
        self.__class__.organism_count += 1
        self._validate_sequence(sequence)
        self.bits, self.rc_bits = self._encode(sequence)
        self.length = len(self)

    def __len__(self):
        return len(self.bits) // 2

    def _validate_sequence(self, sequence: str):
        if not all(base in self._encode_map for base in sequence.upper()):
            bad = sorted({b for b in sequence.upper() if b not in self._encode_map})
            raise ValueError(f"Sequence contains invalid nucleotides: {bad}")

    def _encode(self, sequence: str):
        bits = bitarray()
        rc_bits = bitarray()
        seq = sequence.upper()
        for base in seq:
            bits.extend(self._encode_map[base])
        for base in seq[::-1]:
            rc_bits.extend(self._comp_encode_map[base])
        return bits, rc_bits

    def decode(self) -> str:
        seq = []
        for i in range(self.length):
            two_bits = self.bits[2 * i: 2 * i + 2]
            seq.append(self._decode_map[two_bits.to01()])
        return ''.join(seq)

    def rc_decode(self) -> str:
        seq = []
        for i in range(self.length):
            two_bits = self.rc_bits[2 * i: 2 * i + 2]
            seq.append(self._decode_map[two_bits.to01()])
        return ''.join(seq)

    def get_encoded_slice(self, index: int, slice_length: int, is_rev_comp: bool = False) -> bitarray:
        assert slice_length > 0, "slice_length must be positive"
        assert 0 <= index < self.length, "Index out of bounds"

        if not is_rev_comp:
            if index + slice_length > self.length:
                return bitarray()
            bit_start = 2 * index
            bit_end = 2 * (index + slice_length)
            return self.bits[bit_start:bit_end]
        else:
            if index - slice_length + 1 < 0:
                return bitarray()
            bit_start = 2 * (self.length - index - 1)
            bit_end = 2 * (self.length - index - 1 + slice_length)
            return self.rc_bits[bit_start:bit_end]
