require(Rcpp)
require(RcppEigen)

# Rcpp functions
Rcpp::cppFunction('
Eigen::MatrixXd cpp_mat_vec_prod(const Eigen::MatrixXd X, const Eigen::VectorXd y){
Eigen::VectorXd output;
output = X*y;
return output;
}', depends = 'RcppEigen')

Rcpp::cppFunction('
Eigen::MatrixXd cpp_prod(const Eigen::MatrixXd X, const Eigen::MatrixXd Y){
Eigen::MatrixXd output;
output = X*Y;
return output;
}', depends = 'RcppEigen')

