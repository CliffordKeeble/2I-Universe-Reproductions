# Paper 174 -- icosahedral spectrum: reproduction results

```
Paper 174 -- The Icosahedral Spectrum: m_k reproduction (k = 0..120)
======================================================================
[PASS] Four-method agreement k=0..120 (mismatches: [])
[PASS] Character sum within 1e-10 of integers (max deviation 1.676e-49)
[PASS] Series identity: Klein (1-t^60)/(1-t^30) == Kostant (1+t^30) for all k<=120
[PASS] First non-trivial survivor at k=12 (got k=12)
[PASS] lambda_1 = k(k+2) = 12*14 = 168 (got 168)
[PASS] All odd k killed (alive odd k: [])
[PASS] Exactly 15 killed even modes == expected (got [2, 4, 6, 8, 10, 14, 16, 18, 22, 26, 28, 34, 38, 46, 58])
[PASS] 15 killed even modes, largest = 58 (got count=15, max=58)
[PASS] Frobenius(<6,10,15>) = 29 (got 29)
[PASS] 2 * Frobenius = 58 largest even non-member (got 58)
[PASS] Every even k > 58 survives (dead even k>58: [])
[PASS] m_12 = 1 (got 1)
[PASS] m_20 = 1 (got 1)
[PASS] m_24 = 1 (got 1)
[PASS] m_60 = 2 (got 2)

m_k for k = 0..120 (survivors only shown with their Laplacian eigenvalue k(k+2)):
  k=  0  m_k=1  lambda=k(k+2)=0
  k= 12  m_k=1  lambda=k(k+2)=168
  k= 20  m_k=1  lambda=k(k+2)=440
  k= 24  m_k=1  lambda=k(k+2)=624
  k= 30  m_k=1  lambda=k(k+2)=960
  k= 32  m_k=1  lambda=k(k+2)=1088
  k= 36  m_k=1  lambda=k(k+2)=1368
  k= 40  m_k=1  lambda=k(k+2)=1680
  k= 42  m_k=1  lambda=k(k+2)=1848
  k= 44  m_k=1  lambda=k(k+2)=2024
  k= 48  m_k=1  lambda=k(k+2)=2400
  k= 50  m_k=1  lambda=k(k+2)=2600
  k= 52  m_k=1  lambda=k(k+2)=2808
  k= 54  m_k=1  lambda=k(k+2)=3024
  k= 56  m_k=1  lambda=k(k+2)=3248
  k= 60  m_k=2  lambda=k(k+2)=3720
  k= 62  m_k=1  lambda=k(k+2)=3968
  k= 64  m_k=1  lambda=k(k+2)=4224
  k= 66  m_k=1  lambda=k(k+2)=4488
  k= 68  m_k=1  lambda=k(k+2)=4760
  k= 70  m_k=1  lambda=k(k+2)=5040
  k= 72  m_k=2  lambda=k(k+2)=5328
  k= 74  m_k=1  lambda=k(k+2)=5624
  k= 76  m_k=1  lambda=k(k+2)=5928
  k= 78  m_k=1  lambda=k(k+2)=6240
  k= 80  m_k=2  lambda=k(k+2)=6560
  k= 82  m_k=1  lambda=k(k+2)=6888
  k= 84  m_k=2  lambda=k(k+2)=7224
  k= 86  m_k=1  lambda=k(k+2)=7568
  k= 88  m_k=1  lambda=k(k+2)=7920
  k= 90  m_k=2  lambda=k(k+2)=8280
  k= 92  m_k=2  lambda=k(k+2)=8648
  k= 94  m_k=1  lambda=k(k+2)=9024
  k= 96  m_k=2  lambda=k(k+2)=9408
  k= 98  m_k=1  lambda=k(k+2)=9800
  k=100  m_k=2  lambda=k(k+2)=10200
  k=102  m_k=2  lambda=k(k+2)=10608
  k=104  m_k=2  lambda=k(k+2)=11024
  k=106  m_k=1  lambda=k(k+2)=11448
  k=108  m_k=2  lambda=k(k+2)=11880
  k=110  m_k=2  lambda=k(k+2)=12320
  k=112  m_k=2  lambda=k(k+2)=12768
  k=114  m_k=2  lambda=k(k+2)=13224
  k=116  m_k=2  lambda=k(k+2)=13688
  k=118  m_k=1  lambda=k(k+2)=14160
  k=120  m_k=3  lambda=k(k+2)=14640
Total survivors in 0..120: 46

======================================================================
RESULT: PASS -- all acceptance targets reproduced.
```
