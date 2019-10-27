# BOLTZMANN
boltzmann-machine Many body 

ISING: 이징체인 모델 푼 것

ISING_SR: 이징 모델 풀 때, 스토캐스틱 리컨피겨레이션 이용한 것

Heisenberg 1D: 하이젠버그 체인 문제 푼 것

Heisenberg 2D: 하이젠버그 사각격자 문제 푼 것

ED 붙은 것: Exact diagonalization 코드들. 머신러닝 결과와 비교하기 위해 작성한 코드.

          컴파일 하려면 LAPACK 설치 후 -llapack 옵션 주면 됨.
          
          ISING_ED 코드는 자기장을 입력으로 받게 되어 있음.
          
          Heisenberg는 모두 J=1로 소스안에 설정됨. (스핀크기는 1/2로 되어 있음)
          
* RBM_2SET : 크기와 페이즈 파라메터를 따로 정해서 계산하는 코드. 
   
             Heisenberg 1D 문제의 계산 결과 포함 
             
             계산이 오래 걸리고 대부분 계산이 실패함.
             
* 2D 하이젠버그 격자 제외하고는 10분 이내에 결과가 나옴. ISING의 visible unit = 10, 

                                              Heisenber chain의 visible unit = 5개, 

                                              hidden unit 은 visible unit의 개수와 같게 설정되어 있음.
                                         
 * Heisenber 2D는 3x3까지 해봤는데, 시간이 오래 걸리고 계산이 불안함. 2x2는 계산이 빠르게 나옴.
 
 * Translation symmetry를 이용한 Ising model의 계산코드와, Time evolution code 는 이창우 박사의 서버에 있을 것으로 추정.
