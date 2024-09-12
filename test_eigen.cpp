#include <eigen3/Eigen/Dense>
#include <iostream>

using namespace std;
template <int predH, int controlH> class MPC {
private:
  double kl = 1.0;
  double tl = 1.2;
  double tsample = 0.1;
  double th = 0.8;

  template <typename MatrixType>
  MatrixType pow(const MatrixType &x, int power) {
    MatrixType res;
    for (size_t i = 0; i < res.rows(); i++)
    {
        for (size_t j = 0; j < res.cols(); j++)
        {
            if (i == j)
            { 
                res(i,j) =1;
            }
            else
            {
                res(i,j) =0;
            }
            
        } 
    }
    
    for (size_t i = 0; i < power; i++) {
      res *= x;
    }
    return res;
  }

public:
  MPC(/* args */);
  ~MPC();
};

template <int predH, int controlH> MPC<predH, controlH>::MPC(/* args */) {
  Eigen::Matrix3d A;
  A << 1, tsample, -th * tsample, 0, 1, -tsample, 0, 0, 1 - tsample / tl;

  Eigen::Vector3d B;
  B << 0, 0, tsample * kl / tl;

  Eigen::Vector3d G;
  G << 0, tsample, 0;

  Eigen::Matrix3d C;
  C << 1, 0, 0, 0, 1, 0, 0, 0, 1;

  Eigen::Matrix<Eigen::Matrix3d, predH, 1> Mx;
  for (size_t i = 0; i < predH; i++) {
    Mx(i) = C * this->pow(A, i + 1);
  }
//   cout<<"Mx: "<<endl;
//   for (size_t i = 0; i < predH; i++)
//   {
//     cout<<Mx(i)<<endl<<endl;
//   }
  
  Eigen::Matrix<Eigen::Matrix3d, predH, 1> Me;
  for (size_t i = 0; i < predH; i++) {
    Me(i) = C * this->pow(A, i);
  }

  Eigen::Matrix<Eigen::Vector3d, predH, 1> Mu;
  for (size_t i = 0; i < predH; i++) {
    Eigen::Vector3d t = Eigen::Vector3d::Zero();
    for (size_t j = 0; j < i+1; j++)
    {
        t+=this->pow(A,j)*B;
    }
    Mu(i) = C*t;
  }
//   cout<<"Mu: "<<endl;
//   for (size_t i = 0; i < predH; i++)
//   {
//     cout<<Mu(i)<<endl<<endl;
//   }

  Eigen::Matrix<Eigen::Vector3d, predH, controlH> Mdu;
  for (size_t j = 0; j < controlH; j++) // loop col
  {
    for (size_t k = 0; k < j; k++) // loop row
    {
      Mdu(k, j) = Eigen::Vector3d::Zero();
    }
    for (size_t i = j; i < predH; i++) { // 
      Eigen::Vector3d t = Eigen::Vector3d::Zero();
      for (size_t idx = 0; idx < i - j + 1; idx++) { // 
        t += this->pow(A, idx) * B;
      }
      Mdu(i, j) = C * t;
    }
  }
  cout<<"Mdu: "<<endl;
  for (size_t i = 0; i < predH; i++) {

    cout<<Mdu(i, 0)<<endl<<endl;

  }

  Eigen::Matrix<Eigen::Vector3d, predH, predH + 1> Mv;
  for (size_t j = 0; j < predH + 1; j++) {
    for (size_t k = 0; k < j; k++) // loop row
    {
      Mv(k, j) = Eigen::Vector3d::Zero();
    }
    for (size_t i = j; i < predH; i++) {
      Mv(i, j) = C * this->pow(A, i - j) * G;
    }
  }

}
template <int predH, int controlH> MPC<predH, controlH>::~MPC() {}

int main() {
  Eigen::MatrixXi m(2, 2);
  m << 1, 2, 3, 4;
  cout << m << endl;
  MPC<20,10> mpc;
  return 0;
}