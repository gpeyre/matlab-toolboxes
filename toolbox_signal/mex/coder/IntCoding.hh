// Copyright 1996 John M. Danskin 7/30/96
//
// Permission is granted to use this software for any purpose as
// long as this notice stays attached to this software.
//

class CdeltaEncode {
public:
    CdeltaEncode(BitOut *bo) : output(bo) {};
    void EncodePositive(int i);
    void EncodeNonNegative(int i);
    void Encode(int i);
private: 
    BitOut *output;
};

class CdeltaDecode {
public:
    CdeltaDecode(BitIn *bi) : input(bi) {};
    int Decode(void);
    int DecodePositive(void);
    int DecodeNonNegative(void);
private: 
    BitIn *input;
};
