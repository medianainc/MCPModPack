#include <Rcpp.h>
#include <RcppNumerical.h>

using namespace Numer;
using namespace std;
using namespace Rcpp;

typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

// Global variables
int endpoint_index = 1;
int n_models = 6;
double pi = 3.14159265359;
NumericMatrix temp_matrix;
double final_gradient;
int iter;

// Long vector of overdispersion parameters in the negative binomial distribution (one value for each patient)
vector<double> theta_vector;

// Information required to fit a dose-response model and fitting results
struct ModelInformation 
{
    int flag;
    int n_parameters;
    vector<double> initial_values;
    vector<double> coef;
    int status;
    double criterion;
    double target_dose;
    double convergence_criterion;

};

// Information on estimated model parameters in the simulation module
struct ModelParameters 
{
    int scenario;
    int model;
    vector<double> coef;
};

// Vector of optimal contrast matrices
struct OptimalContrastMatrix 
{
    NumericMatrix matrix;
};

// Vector of constant values
vector<double> FillVec(const int &n, const double &value) {

    vector<double> result(n); 
    for (int i = 0; i < n; i++) {
        result[i] = value;
    }
    return result;  
}

// Quantile of the gamma distribution
double rcpp_qgamma(const double &x, const double &shape, const double &scale) {
    NumericVector vec_input(1), vec_output(1);
    vec_input[0] = x;
    vec_output = Rcpp::qgamma(vec_input, shape, scale);
    return vec_output[0];
}

// This function is based on the code written by John Cook, see https://www.johndcook.com/blog/csharp_log_factorial/
double LogFactorial(const double &n)
{

    if (n > 254)
    {
        double x = n + 1;
        return (x - 0.5)*log(x) - x + 0.5*log(2*pi) + 1.0/(12.0*x);
    }
    else
    {
        double lf[] =        
        {
            0.000000000000000,
            0.000000000000000,
            0.693147180559945,
            1.791759469228055,
            3.178053830347946,
            4.787491742782046,
            6.579251212010101,
            8.525161361065415,
            10.604602902745251,
            12.801827480081469,
            15.104412573075516,
            17.502307845873887,
            19.987214495661885,
            22.552163853123421,
            25.191221182738683,
            27.899271383840894,
            30.671860106080675,
            33.505073450136891,
            36.395445208033053,
            39.339884187199495,
            42.335616460753485,
            45.380138898476908,
            48.471181351835227,
            51.606675567764377,
            54.784729398112319,
            58.003605222980518,
            61.261701761002001,
            64.557538627006323,
            67.889743137181526,
            71.257038967168000,
            74.658236348830158,
            78.092223553315307,
            81.557959456115029,
            85.054467017581516,
            88.580827542197682,
            92.136175603687079,
            95.719694542143202,
            99.330612454787428,
            102.968198614513810,
            106.631760260643450,
            110.320639714757390,
            114.034211781461690,
            117.771881399745060,
            121.533081515438640,
            125.317271149356880,
            129.123933639127240,
            132.952575035616290,
            136.802722637326350,
            140.673923648234250,
            144.565743946344900,
            148.477766951773020,
            152.409592584497350,
            156.360836303078800,
            160.331128216630930,
            164.320112263195170,
            168.327445448427650,
            172.352797139162820,
            176.395848406997370,
            180.456291417543780,
            184.533828861449510,
            188.628173423671600,
            192.739047287844900,
            196.866181672889980,
            201.009316399281570,
            205.168199482641200,
            209.342586752536820,
            213.532241494563270,
            217.736934113954250,
            221.956441819130360,
            226.190548323727570,
            230.439043565776930,
            234.701723442818260,
            238.978389561834350,
            243.268849002982730,
            247.572914096186910,
            251.890402209723190,
            256.221135550009480,
            260.564940971863220,
            264.921649798552780,
            269.291097651019810,
            273.673124285693690,
            278.067573440366120,
            282.474292687630400,
            286.893133295426990,
            291.323950094270290,
            295.766601350760600,
            300.220948647014100,
            304.686856765668720,
            309.164193580146900,
            313.652829949878990,
            318.152639620209300,
            322.663499126726210,
            327.185287703775200,
            331.717887196928470,
            336.261181979198450,
            340.815058870798960,
            345.379407062266860,
            349.954118040770250,
            354.539085519440790,
            359.134205369575340,
            363.739375555563470,
            368.354496072404690,
            372.979468885689020,
            377.614197873918670,
            382.258588773060010,
            386.912549123217560,
            391.575988217329610,
            396.248817051791490,
            400.930948278915760,
            405.622296161144900,
            410.322776526937280,
            415.032306728249580,
            419.750805599544780,
            424.478193418257090,
            429.214391866651570,
            433.959323995014870,
            438.712914186121170,
            443.475088120918940,
            448.245772745384610,
            453.024896238496130,
            457.812387981278110,
            462.608178526874890,
            467.412199571608080,
            472.224383926980520,
            477.044665492585580,
            481.872979229887900,
            486.709261136839360,
            491.553448223298010,
            496.405478487217580,
            501.265290891579240,
            506.132825342034830,
            511.008022665236070,
            515.890824587822520,
            520.781173716044240,
            525.679013515995050,
            530.584288294433580,
            535.496943180169520,
            540.416924105997740,
            545.344177791154950,
            550.278651724285620,
            555.220294146894960,
            560.169054037273100,
            565.124881094874350,
            570.087725725134190,
            575.057539024710200,
            580.034272767130800,
            585.017879388839220,
            590.008311975617860,
            595.005524249382010,
            600.009470555327430,
            605.020105849423770,
            610.037385686238740,
            615.061266207084940,
            620.091704128477430,
            625.128656730891070,
            630.172081847810200,
            635.221937855059760,
            640.278183660408100,
            645.340778693435030,
            650.409682895655240,
            655.484856710889060,
            660.566261075873510,
            665.653857411105950,
            670.747607611912710,
            675.847474039736880,
            680.953419513637530,
            686.065407301994010,
            691.183401114410800,
            696.307365093814040,
            701.437263808737160,
            706.573062245787470,
            711.714725802289990,
            716.862220279103440,
            722.015511873601330,
            727.174567172815840,
            732.339353146739310,
            737.509837141777440,
            742.685986874351220,
            747.867770424643370,
            753.055156230484160,
            758.248113081374300,
            763.446610112640200,
            768.650616799717000,
            773.860102952558460,
            779.075038710167410,
            784.295394535245690,
            789.521141208958970,
            794.752249825813460,
            799.988691788643450,
            805.230438803703120,
            810.477462875863580,
            815.729736303910160,
            820.987231675937890,
            826.249921864842800,
            831.517780023906310,
            836.790779582469900,
            842.068894241700490,
            847.352097970438420,
            852.640365001133090,
            857.933669825857460,
            863.231987192405430,
            868.535292100464630,
            873.843559797865740,
            879.156765776907600,
            884.474885770751830,
            889.797895749890240,
            895.125771918679900,
            900.458490711945270,
            905.796028791646340,
            911.138363043611210,
            916.485470574328820,
            921.837328707804890,
            927.193914982476710,
            932.555207148186240,
            937.921183163208070,
            943.291821191335660,
            948.667099599019820,
            954.046996952560450,
            959.431492015349480,
            964.820563745165940,
            970.214191291518320,
            975.612353993036210,
            981.015031374908400,
            986.422203146368590,
            991.833849198223450,
            997.249949600427840,
            1002.670484599700300,
            1008.095434617181700,
            1013.524780246136200,
            1018.958502249690200,
            1024.396581558613400,
            1029.838999269135500,
            1035.285736640801600,
            1040.736775094367400,
            1046.192096209724900,
            1051.651681723869200,
            1057.115513528895000,
            1062.583573670030100,
            1068.055844343701400,
            1073.532307895632800,
            1079.012946818975000,
            1084.497743752465600,
            1089.986681478622400,
            1095.479742921962700,
            1100.976911147256000,
            1106.478169357800900,
            1111.983500893733000,
            1117.492889230361000,
            1123.006317976526100,
            1128.523770872990800,
            1134.045231790853000,
            1139.570684729984800,
            1145.100113817496100,
            1150.633503306223700,
            1156.170837573242400,
        };
        return lf[(int)n];
    }
}

// Square a number
double Sq(const double &x) {
    return x * x ;
}

// Truncated natural log
double TruncLog(const double &x) {
    double trunc_log = log(0.001);
    if (x > 0.001) trunc_log = log(x);
    return trunc_log;
}

double Logit(const double &x) {
    return log(x /(1.0 - x));
}

double AntiLogit(const double &x) {
    double result;
    if (x <= 0.0) result = exp(x) /(1.0 + exp(x)); else result = 1.0 /(1.0 + exp(-x));
    return result;
}

// Vector of binary values 
vector<double> rcpp_binary(const int &n, const double &prop) {

    NumericVector result = Rcpp::rbinom(n, 1, prop);
    vector<double> result_vector = as<vector<double>>(result);
    return result_vector;  

}

// Vector of uniformly distributed values 
vector<double> rcpp_uniform(const int &n, const double &min, const double &max) {

    NumericVector result = Rcpp::runif(n, min, max);
    vector<double> result_vector = as<vector<double>>(result);
    return result_vector;  

}

// Vector of normally distributed values 
vector<double> rcpp_normal(const int &n, const double &mean, const double &sd) {

    NumericVector result = Rcpp::rnorm(n, mean, sd);
    vector<double> result_vector = as<vector<double>>(result);
    return result_vector;  

}

// Vector of negative binomial distributed values 
vector<double> rcpp_nbinom(const int &n, const double &r, const double &p) {

    NumericVector result = Rcpp::rnbinom(n, r, p);
    vector<double> result_vector = as<vector<double>>(result);
    return result_vector;  

}

// Compute averages by dividing by the number of simulations 
vector<double> AverageSimResultsVector(const vector<double> &vec, const int &nsims) {

    int i;
    int m = vec.size();
    vector<double> ave(m); 

    for (i = 0; i < m; i++) {
        ave[i] = vec[i] / nsims;
    }

    return ave;  

}

double Quantile(const vector<double> &vec, const double &prob) {

    int i, m = vec.size();

    Rcpp::NumericVector rvec(m);

    for(i = 0; i < m; i++) rvec[i] = vec[i];

    // Obtain environment containing function
    Rcpp::Environment base("package:stats"); 

    // Make function callable from C++
    Rcpp::Function rquantile = base["quantile"];    

    // Call the function and receive its list output
    Rcpp::NumericVector res = rquantile(Rcpp::_["x"] = rvec,
                                      Rcpp::_["probs"] = prob,
                                      Rcpp::_["na.rm"] = true
                                      ); 

    return res[0];

}

vector<double> ConvertToDoubleVector(const NumericVector &numeric_vec) {

    int i, m = numeric_vec.size();
    vector<double> vec(m);

    for(i = 0; i < m; i++) vec[i] = numeric_vec[i];

    return vec;    

}

// Compute averages by dividing by the number of simulations 
NumericMatrix AverageSimResultsMatrix(const NumericMatrix &mat, const int &nsims) {

    int i, j;
    int m1 = mat.nrow();
    int m2 = mat.ncol();
    NumericMatrix ave(m1, m2); 

    for (i = 0; i < m1; i++) {
        for (j = 0; j < m2; j++) {
            ave(i, j) = mat(i, j) / nsims;
        }
    }

    return ave;  

}

// Compute the mean of a vector
double MeanVec(const vector<double> &vec)
{

   return(accumulate(vec.begin(), vec.end(), 0.0) / vec.size());

}

double MaxVec(const vector<double> &x) {
    return *std::max_element(x.begin(), x.end());
}

vector<double> AddVec(const vector<double> &x, const vector<double> &y) {
    int i, m = x.size();
    vector<double> sum(m);
    for(i = 0; i < m; ++i) sum[i] = x[i] + y[i];
    return sum;
}

vector<double> FitLinearModel(const vector<double> &x, const vector<double> &y) {

    int i, n = x.size();
    vector<double> coef(2);
    double temp1 = 0.0, temp2 = 0.0;

    double x_mean = MeanVec(x);
    double y_mean = MeanVec(y);

    for(i = 0; i < n; i++) {
        temp1 += (x[i] - x_mean) * (y[i] - y_mean);
        temp2 += Sq(x[i] - x_mean);
    }

    coef[1] = temp1 / temp2;
    coef[0] = y_mean - coef[1] * x_mean;

    return coef;  
}


double DoseResponseFunction(const double &x, const int &model, const vector<double> &beta, const double &direction_index) {

    double y = 0.0, den;

    // Linear model
    if (model == 1) {
        y = beta[0] + beta[1] * x;
    }

    // Quadratic model
    if (model == 2) {
        y = beta[0] + beta[1] * x + beta[2] * Sq(x);
    }

    // Exponential model
    if (model == 3) {
        y = beta[0] + beta[1] * (exp(x / beta[2]) - 1.0);
    }

    // Emax model
    if (model == 4) {
        y = beta[0] + beta[1] * x / (beta[2] + x);
    }

    // Logistic model
    if (model == 5) {
        den = 1.0 + exp((beta[2] - x) / beta[3]);
        y = beta[0] + beta[1] / den;
    }

    // SigEmax model
    if (model == 6) {
        // temp = max(beta[2], 0.0);
        den = pow(x, beta[3]) + pow(beta[2], beta[3]);
        y = beta[0] + beta[1] * pow(x, beta[3]) / den;
    }

    // Binary endpoint
    if (endpoint_index == 2) {
        y = AntiLogit(y);
    }

    // Count endpoint
    if (endpoint_index == 3) {
        y = exp(y);
    }

    return y;


}

// Compute the model parameters to match the placebo and maximum effects
vector<double> ComputeDoseResponseFunctionParameters(const int &model, const double &placebo_effect, const double &max_effect, const double &max_dose, const vector<double> &nonlinear_parameters) {

    vector<double> coef(4), temp_coef(4);
    double temp = 0.0;

    // Quadratic model
    if (model == 2) {
        temp = - 1.0 / (2.0 * nonlinear_parameters[0]);
        coef[0] = placebo_effect;
        coef[1] = 2.0 * max_effect / max_dose;
        coef[2] = - 0.5 * coef[1] / max_dose;
    }

    // Exponential model
    if (model == 3) {
        coef[0] = placebo_effect;
        coef[1] = max_effect / (exp(max_dose / nonlinear_parameters[0]) - 1.0);
        coef[2] = nonlinear_parameters[0];
    }

    // Emax model
    if (model == 4) {
        coef[0] = placebo_effect;
        coef[1] = max_effect * (nonlinear_parameters[0] + max_dose) / max_dose;
        coef[2] = nonlinear_parameters[0];
    }

    // Logistic model
    if (model == 5) {
        temp_coef[0] = 0.0;
        temp_coef[1] = 1.0;
        temp_coef[2] = nonlinear_parameters[0];
        temp_coef[3] = nonlinear_parameters[1];
        temp =  max_effect / (DoseResponseFunction(max_dose, 5, temp_coef, 1) - DoseResponseFunction(0.0, 5, temp_coef, 1));
        coef[0] = placebo_effect- temp * DoseResponseFunction(0.0, 5, temp_coef, 1);
        coef[1] = temp;
        coef[2] = nonlinear_parameters[0];
        coef[3] = nonlinear_parameters[1];
    }

    // SigEmax model
    if (model == 6) {
        coef[0] = placebo_effect;
        coef[1] = max_effect * (pow(nonlinear_parameters[0], nonlinear_parameters[1]) + pow(max_dose, nonlinear_parameters[1])) / pow(max_dose, nonlinear_parameters[1]);
        coef[2] = nonlinear_parameters[0];
        coef[3] = nonlinear_parameters[1];
    }

    return coef;

}


double FindMaxEffect(const int &model, const vector<double> &beta, const int &max_dose, const double &direction_index) {

    double max_effect = 0.0, dose;

    // All models but quadratic model
    if (model != 2) {

        max_effect = DoseResponseFunction(max_dose, model, beta, direction_index); 

    } else {

        // Quadratic model: Find the dose corresponding to the max effect
        dose = - beta[1] / (2.0 * beta[2]);

        if (dose < 0.0) dose = 0.0;
        if (dose > max_dose) dose = max_dose;

        max_effect = DoseResponseFunction(dose, model, beta, direction_index); 

    }

    return max_effect;


}

// Compute the target dose
double FindTargetDose(const int &model, const vector<double> &beta, const double &delta, const double &direction_index) {

  double dose = 0.0, local_delta, d1, d2, d3;

  local_delta = delta;

  // Compute the delta on the endpoint-specific scale
  if (endpoint_index == 1) local_delta = delta;

  if (endpoint_index == 2) {

    if (DoseResponseFunction(0.0, model, beta, 1.0) + delta <= 0.0 || DoseResponseFunction(0.0, model, beta, direction_index) + delta >= 1.0) {
        dose = -1.0;
        return dose;
    } else {
        local_delta = Logit(DoseResponseFunction(0.0, model, beta, 1.0) + delta) - Logit(DoseResponseFunction(0.0, model, beta, 1.0));
    }

  }

  if (endpoint_index == 3) {

    if (DoseResponseFunction(0.0, model, beta, 1.0) + delta <= 0.0) {
        dose = -1.0;
        return dose;
    } else {
        local_delta = log(DoseResponseFunction(0.0, model, beta, 1.0) + delta) - log(DoseResponseFunction(0.0, model, beta, 1.0));
    }

  }


  if (dose >= 0.0) {

      // Linear
      if (model == 1) {

        if (abs(beta[1]) <= 0.0001) {
          dose = -1.0; 
        } else {
          dose = local_delta / beta[1];
        }
      }

      // Quadratic model
      if (model == 2) {
          if (4.0 * beta[2] * local_delta + Sq(beta[1]) < 0) {
              dose = -1.0;
          } else {
              d1 = -(sqrt(4.0 * beta[2] * local_delta + Sq(beta[1])) + beta[1])/(2.0 * beta[2]);
              d2 = (sqrt(4.0 * beta[2] * local_delta + Sq(beta[1])) - beta[1])/(2.0 * beta[2]);
              if (d1 <= 0.0 && d2 <= 0.0) dose = -1.0;            
              if (d1 <= 0.0 && d2 > 0.0) dose = d2;            
              if (d1 > 0.0 && d2 <= 0.0) dose = d1;            
              if (d1 > 0.0 && d2 > 0.0) dose = min(d1, d2);            

          }
      }

      // Exponential
      if (model == 3) {
          if (direction_index == 1) {
            if (beta[1] > 0 &&  local_delta + beta[1] > 0.0) {
                dose = beta[2] * (-log(beta[1]) + log(local_delta + beta[1]));  
            } else {
                dose = -1.0;   
            } 
          }

          if (direction_index == -1) {
            if (beta[1] < 0 && local_delta + beta[1] < 0.0) {
                dose = beta[2] * (-log(-beta[1]) + log(-local_delta - beta[1]));  
            } else {
                dose = -1.0;   
            } 
          }

      }

      // Emax
      if (model == 4) {
        if (direction_index == 1 && local_delta >= beta[1]) {      
             dose = -1.0; 
        }
        if (direction_index == -1 && local_delta <= beta[1]) {      
             dose = -1.0; 
        }
        if (dose == 0) {      
            dose = local_delta * beta[2]/(beta[1]-local_delta);
        }
      }


    // Logistic model
    if (model == 5) {
      d1 = 1.0 /(1.0 + exp(beta[2]/beta[3])); 
      if (direction_index == 1 && local_delta >= beta[1] * (1.0 - d1)) {
            dose = -1.0; 
        } 
      if (direction_index == -1 && local_delta <= beta[1] * (1.0 - d1)) {
            dose = -1.0; 
        } 

      if (dose == 0) {
            d1 = exp(beta[2]/beta[3]);
            d2 = d1 * beta[1] - local_delta * d1 - local_delta;
            d3 = beta[1] + local_delta * d1 + local_delta;
            dose = beta[2] - beta[3] * log(d2 / d3);
        }
    }

    // SigEmax model
    if (model == 6) {
      if (direction_index == 1 && local_delta >= beta[1]) {              
        dose = -1.0; 
      }
      if (direction_index == -1 && local_delta <= beta[1]) {
        dose = -1.0; 
      }
      if (dose == 0) {      
        dose = pow(local_delta * pow(beta[2], beta[3])/(beta[1] - local_delta), 1.0 / beta[3]);
      }

    }
  }

  if (dose <= 0.0) dose = -1.0;
  if (dose >= 10000.0) dose = -1.0;

  return dose;

}


class RegressionLinear: public MFuncGrad
{
private:
    const MapVec X;
    const MapVec Y;
public:
    RegressionLinear(const MapVec x_, const MapVec y_) : X(x_), Y(y_) {}

    double f_grad(Constvec& beta, Refvec grad)
    {

        int i, n = Y.size(), n_parameters;
        double neg_log_likelihood = 0.0, p, mu, e0, delta, sigma, main_term;
        Eigen::VectorXd gradient; 

        // Normal endpoint
        if (endpoint_index == 1) {

            n_parameters = 3;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;           

            // Restricted estimates
            e0 = beta[0];    
            delta = beta[1];    
            sigma = max(beta[2], 0.0001);    

            for(i = 0; i < n; i++) {
                mu = e0 + delta * X[i];
                neg_log_likelihood += (log(sqrt(2 * pi) * sigma) + Sq(Y[i] - mu) / (2 * Sq(sigma)));
                main_term = (mu - Y[i]) / Sq(sigma);
                gradient[0] += main_term;
                gradient[1] += main_term * X[i];
                gradient[2] += (1.0 / sigma - Sq(Y[i] - mu) / (sigma * sigma * sigma));
            }

        }

        // Binary endpoint
        if (endpoint_index == 2) {

            n_parameters = 2;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;           

            // Restricted estimates
            e0 = beta[0];    
            delta = beta[1];    

            for(i = 0; i < n; i++) {
                mu = e0 + delta * X[i];

                if (mu <= 0.0) {
                    p = exp(mu) / (1.0 + exp(mu));
                    neg_log_likelihood += (log(1.0 + exp(mu)) - Y[i] * mu);
                } else {
                    p = 1.0 / (1.0 + exp(-mu));
                    neg_log_likelihood += (mu + log(1.0 + exp(-mu)) - Y[i] * mu);
                } 

                main_term = p - Y[i];
                gradient[0] += main_term;
                gradient[1] += main_term * X[i];
            }

        }

        // Count endpoint
        if (endpoint_index == 3) {

            n_parameters = 2;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;           

            // Restricted estimates
            e0 = beta[0];    
            delta = beta[1];    

            for(i = 0; i < n; i++) {
                mu = e0 + delta * X[i];
                if (mu <= 0.0) {
                    neg_log_likelihood += (lgamma(theta_vector[i]) + LogFactorial(Y[i]) - lgamma(theta_vector[i] + Y[i]) + (theta_vector[i] + Y[i]) * log(theta_vector[i] + exp(mu)) - Y[i] * mu - theta_vector[i] * log(theta_vector[i]));
                    main_term = theta_vector[i] * (exp(mu) - Y[i]) / (theta_vector[i] + exp(mu));
                } else {
                    neg_log_likelihood += (lgamma(theta_vector[i]) + LogFactorial(Y[i]) - lgamma(theta_vector[i] + Y[i]) + (theta_vector[i] + Y[i]) * (mu + log(theta_vector[i] * exp(-mu) + 1.0)) - Y[i] * mu - theta_vector[i] * log(theta_vector[i]));
                    main_term = theta_vector[i] * (1.0 - Y[i] * exp(-mu)) / (theta_vector[i] * exp(-mu) + 1.0);                    
                }
                gradient[0] += main_term;
                gradient[1] += main_term * X[i];
            }

        }        

        final_gradient = 0.0;
        for (i = 0; i < n_parameters; i++) final_gradient += abs(gradient[i]);

        const double f = neg_log_likelihood;
        grad.noalias() = gradient;

        return f;
    }
};


class RegressionQuadratic: public MFuncGrad
{
private:
    const MapVec X;
    const MapVec Y;
public:
    RegressionQuadratic(const MapVec x_, const MapVec y_) : X(x_), Y(y_) {}

    double f_grad(Constvec& beta, Refvec grad)
    {

        int i, n = Y.size(), n_parameters;
        double neg_log_likelihood = 0.0, p, mu, e0, delta1, delta2, sigma, main_term;
        Eigen::VectorXd gradient; 

        // Normal endpoint
        if (endpoint_index == 1) {

            n_parameters = 4;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;           

            // Restricted estimates
            e0 = beta[0];    
            delta1 = beta[1];    
            delta2 = beta[2];    
            sigma = max(beta[3], 0.0001);    

            for(i = 0; i < n; i++) {
                mu = e0 + delta1 * X[i] + delta2 * Sq(X[i]);
                neg_log_likelihood += (log(sqrt(2 * pi) * sigma) + Sq(Y[i] - mu) / (2 * Sq(sigma)));
                main_term = (mu - Y[i]) / Sq(sigma);
                gradient[0] += main_term;
                gradient[1] += main_term * X[i];
                gradient[2] += main_term * Sq(X[i]);
                gradient[3] += (1.0 / sigma - Sq(Y[i] - mu) / (sigma * sigma * sigma));
            }

        }

        // Binary endpoint
        if (endpoint_index == 2) {

            n_parameters = 3;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;           

            // Restricted estimates
            e0 = beta[0];    
            delta1 = beta[1];    
            delta2 = beta[2];    

            for(i = 0; i < n; i++) {
                mu = e0 + delta1 * X[i] + delta2 * Sq(X[i]);
                if (mu <= 0.0) {
                    p = exp(mu) / (1.0 + exp(mu));
                    neg_log_likelihood += (log(1.0 + exp(mu)) - Y[i] * mu);
                } else {
                    p = 1.0 / (1.0 + exp(-mu));
                    neg_log_likelihood += (mu + log(1.0 + exp(-mu)) - Y[i] * mu);
                } 
                main_term = p - Y[i];
                gradient[0] += main_term;
                gradient[1] += main_term * X[i];
                gradient[2] += main_term * Sq(X[i]);
            }

        }

        // Count endpoint
        if (endpoint_index == 3) {

            n_parameters = 3;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;           

            // Restricted estimates
            e0 = beta[0];    
            delta1 = beta[1];    
            delta2 = beta[2];    

            for(i = 0; i < n; i++) {
                mu = e0 + delta1 * X[i] + delta2 * Sq(X[i]);
                if (mu <= 0.0) {
                    neg_log_likelihood += (lgamma(theta_vector[i]) + LogFactorial(Y[i]) - lgamma(theta_vector[i] + Y[i]) + (theta_vector[i] + Y[i]) * log(theta_vector[i] + exp(mu)) - Y[i] * mu - theta_vector[i] * log(theta_vector[i]));
                    main_term = theta_vector[i] * (exp(mu) - Y[i]) / (theta_vector[i] + exp(mu));
                } else {
                    neg_log_likelihood += (lgamma(theta_vector[i]) + LogFactorial(Y[i]) - lgamma(theta_vector[i] + Y[i]) + (theta_vector[i] + Y[i]) * (mu + log(theta_vector[i] * exp(-mu) + 1.0)) - Y[i] * mu - theta_vector[i] * log(theta_vector[i]));
                    main_term = theta_vector[i] * (1.0 - Y[i] * exp(-mu)) / (theta_vector[i] * exp(-mu) + 1.0);                    
                }
                gradient[0] += main_term;
                gradient[1] += main_term * X[i];
                gradient[2] += main_term * Sq(X[i]);

            }

        }

        final_gradient = 0.0;
        for (i = 0; i < n_parameters; i++) final_gradient += abs(gradient[i]);

        const double f = neg_log_likelihood;
        grad.noalias() = gradient;

        return f;
    }
};

class RegressionExponential: public MFuncGrad
{
private:
    const MapVec X;
    const MapVec Y;
public:
    RegressionExponential(const MapVec x_, const MapVec y_) : X(x_), Y(y_) {}

    double f_grad(Constvec& beta, Refvec grad)
    {

        int i, n = Y.size(), n_parameters;
        double neg_log_likelihood = 0.0, p, mu, e0, e1, delta, sigma, main_term;
        Eigen::VectorXd gradient; 

        // Normal endpoint
        if (endpoint_index == 1) {

            n_parameters = 4;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;   

            // Restricted estimates
            e0 = beta[0];    
            e1 = beta[1];    
            delta = max(beta[2], 0.01);    
            sigma = max(beta[3], 0.0001);    

            for(i = 0; i < n; i++) {
                mu = e0 + e1 * (exp(X[i] / delta) - 1.0);
                neg_log_likelihood += (log(sqrt(2 * pi) * sigma) + Sq(Y[i] - mu) / (2 * Sq(sigma)));
                main_term = (mu - Y[i]) / Sq(sigma);
                gradient[0] += main_term;
                gradient[1] += main_term * (exp(X[i] / delta) - 1.0);
                gradient[2] += main_term * (- e1 * exp(X[i] / delta) / Sq(delta));
                gradient[3] += (1.0 / sigma - Sq(Y[i] - mu) / (sigma * sigma * sigma));
            }
        }

        // Binary endpoint
        if (endpoint_index == 2) {

            n_parameters = 3;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;           

            // Restricted estimates
            e0 = beta[0];    
            e1 = beta[1];    
            delta = max(beta[2], 0.01);    

            for(i = 0; i < n; i++) {
                mu = e0 + e1 * (exp(X[i] / delta) - 1.0);
                if (mu <= 0.0) {
                    p = exp(mu) / (1.0 + exp(mu));
                    neg_log_likelihood += (log(1.0 + exp(mu)) - Y[i] * mu);
                } else {
                    p = 1.0 / (1.0 + exp(-mu));
                    neg_log_likelihood += (mu + log(1.0 + exp(-mu)) - Y[i] * mu);
                } 
                main_term = p - Y[i];
                gradient[0] += main_term;
                gradient[1] += main_term * (exp(X[i] / delta) - 1.0);
                gradient[2] += main_term * (- e1 * exp(X[i] / delta) / Sq(delta));
            }

        }

        // Count endpoint
        if (endpoint_index == 3) {

            n_parameters = 3;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;   

            // Restricted estimates
            e0 = beta[0];    
            e1 = beta[1];    
            delta = max(beta[2], 0.01);    

            for(i = 0; i < n; i++) {
                mu = e0 + e1 * (exp(X[i] / delta) - 1.0);
                if (mu <= 0.0) {
                    neg_log_likelihood += (lgamma(theta_vector[i]) + LogFactorial(Y[i]) - lgamma(theta_vector[i] + Y[i]) + (theta_vector[i] + Y[i]) * log(theta_vector[i] + exp(mu)) - Y[i] * mu - theta_vector[i] * log(theta_vector[i]));
                    main_term = theta_vector[i] * (exp(mu) - Y[i]) / (theta_vector[i] + exp(mu));
                } else {
                    neg_log_likelihood += (lgamma(theta_vector[i]) + LogFactorial(Y[i]) - lgamma(theta_vector[i] + Y[i]) + (theta_vector[i] + Y[i]) * (mu + log(theta_vector[i] * exp(-mu) + 1.0)) - Y[i] * mu - theta_vector[i] * log(theta_vector[i]));
                    main_term = theta_vector[i] * (1.0 - Y[i] * exp(-mu)) / (theta_vector[i] * exp(-mu) + 1.0);                    
                }
                gradient[0] += main_term;
                gradient[1] += main_term * (exp(X[i] / delta) - 1.0);
                gradient[2] += main_term * (- e1 * exp(X[i] / delta) / Sq(delta));

            }
        }

        final_gradient = 0.0;
        for (i = 0; i < n_parameters; i++) final_gradient += abs(gradient[i]);

        const double f = neg_log_likelihood;
        grad.noalias() = gradient;

        return f;
    }
};

class RegressionEmax: public MFuncGrad
{
private:
    const MapVec X;
    const MapVec Y;
public:
    RegressionEmax(const MapVec x_, const MapVec y_) : X(x_), Y(y_) {}

    double f_grad(Constvec& beta, Refvec grad)
    {

        int i, n = Y.size(), n_parameters;
        double neg_log_likelihood = 0.0, p, mu, e0, emax, ed50, sigma, main_term;
        Eigen::VectorXd gradient;

        // Normal endpoint
        if (endpoint_index == 1) {

            n_parameters = 4;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;   

            // Restricted estimates
            e0 = beta[0];    
            emax = beta[1];    
            ed50 = max(beta[2], 0.01);    
            sigma = max(beta[3], 0.0001);    

            for(i = 0; i < n; i++) {
                mu = e0 + emax * X[i] / (ed50 + X[i]);
                neg_log_likelihood += (log(sqrt(2 * pi) * sigma) + Sq(Y[i] - mu) / (2 * Sq(sigma)));
                main_term = (mu - Y[i]) / Sq(sigma);
                gradient[0] += main_term;
                gradient[1] += main_term * X[i] / (ed50 + X[i]);
                gradient[2] += main_term * (- emax * X[i] / Sq(ed50 + X[i]));
                gradient[3] += (1.0 / sigma - Sq(Y[i] - mu) / (sigma * sigma * sigma));
            }

        }

        // Binary endpoint
        if (endpoint_index == 2) {

            n_parameters = 3;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;           

            // Restricted estimates
            e0 = beta[0];    
            emax = beta[1];    
            ed50 = max(beta[2], 0.01);    

            for(i = 0; i < n; i++) {
                mu = e0 + emax * X[i] / (ed50 + X[i]);              
                if (mu <= 0.0) {
                    p = exp(mu) / (1.0 + exp(mu));
                    neg_log_likelihood += (log(1.0 + exp(mu)) - Y[i] * mu);
                } else {
                    p = 1.0 / (1.0 + exp(-mu));
                    neg_log_likelihood += (mu + log(1.0 + exp(-mu)) - Y[i] * mu);
                } 
                main_term = p - Y[i];
                gradient[0] += main_term;
                gradient[1] += main_term * X[i] / (ed50 + X[i]);
                gradient[2] += main_term * (- emax * X[i] / Sq(ed50 + X[i]));
            }

        }

        // Count endpoint
        if (endpoint_index == 3) {

            n_parameters = 3;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;   

            // Restricted estimates
            e0 = beta[0];    
            emax = beta[1];    
            ed50 = max(beta[2], 0.01);    

            for(i = 0; i < n; i++) {
                mu = e0 + emax * X[i] / (ed50 + X[i]);
                if (mu <= 0.0) {
                    neg_log_likelihood += (lgamma(theta_vector[i]) + LogFactorial(Y[i]) - lgamma(theta_vector[i] + Y[i]) + (theta_vector[i] + Y[i]) * log(theta_vector[i] + exp(mu)) - Y[i] * mu - theta_vector[i] * log(theta_vector[i]));
                    main_term = theta_vector[i] * (exp(mu) - Y[i]) / (theta_vector[i] + exp(mu));
                } else {
                    neg_log_likelihood += (lgamma(theta_vector[i]) + LogFactorial(Y[i]) - lgamma(theta_vector[i] + Y[i]) + (theta_vector[i] + Y[i]) * (mu + log(theta_vector[i] * exp(-mu) + 1.0)) - Y[i] * mu - theta_vector[i] * log(theta_vector[i]));
                    main_term = theta_vector[i] * (1.0 - Y[i] * exp(-mu)) / (theta_vector[i] * exp(-mu) + 1.0);                    
                }
                gradient[0] += main_term;
                gradient[1] += main_term * X[i] / (ed50 + X[i]);
                gradient[2] += main_term * (- emax * X[i] / Sq(ed50 + X[i]));
            }

        }

        final_gradient = 0.0;
        for (i = 0; i < n_parameters; i++) final_gradient += abs(gradient[i]);

        const double f = neg_log_likelihood;
        grad.noalias() = gradient;

        return f;
    }
};

class RegressionLogistic: public MFuncGrad
{
private:
    const MapVec X;
    const MapVec Y;
public:
    RegressionLogistic(const MapVec x_, const MapVec y_) : X(x_), Y(y_) {}

    double f_grad(Constvec& beta, Refvec grad)
    {

        int i, n = Y.size(), n_parameters;
        double neg_log_likelihood = 0.0, p, mu, e0, emax, ed50, delta, sigma, den, main_term;
        Eigen::VectorXd gradient; 

        // Normal endpoint
        if (endpoint_index == 1) {

            n_parameters = 5;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;   

            // Restricted estimates
            e0 = beta[0];    
            emax = beta[1];    
            ed50 = beta[2];    
            delta = beta[3];    
            sigma = max(beta[4], 0.0001);    

            for(i = 0; i < n; i++) {
                den = 1.0 + exp((ed50 - X[i]) / delta);
                mu = e0 + emax / den;
                neg_log_likelihood += (log(sqrt(2 * pi) * sigma) + Sq(Y[i] - mu) / (2 * Sq(sigma)));
                main_term = (mu - Y[i]) / Sq(sigma);
                gradient[0] += main_term;
                gradient[1] += main_term / den;
                gradient[2] += main_term * (- emax * (den - 1.0) / (delta * Sq(den)));
                gradient[3] += main_term * (emax * (den - 1.0) * (ed50 - X[i]) / (Sq(delta * den)));
                gradient[4] += (1.0 / sigma - Sq(Y[i] - mu) / (sigma * sigma * sigma));
            }
        }

        // Binary endpoint
        if (endpoint_index == 2) {

            n_parameters = 4;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;           

            // Restricted estimates
            e0 = beta[0];    
            emax = beta[1];    
            ed50 = beta[2];    
            delta = beta[3];    

            for(i = 0; i < n; i++) {
                den = 1.0 + exp((ed50 - X[i]) / delta);
                mu = e0 + emax / den;
                if (mu <= 0.0) {
                    p = exp(mu) / (1.0 + exp(mu));
                    neg_log_likelihood += (log(1.0 + exp(mu)) - Y[i] * mu);
                } else {
                    p = 1.0 / (1.0 + exp(-mu));
                    neg_log_likelihood += (mu + log(1.0 + exp(-mu)) - Y[i] * mu);
                } 
                main_term = p - Y[i];
                gradient[0] += main_term;
                gradient[1] += main_term / den;
                gradient[2] += main_term * (- emax * (den - 1.0) / (delta * Sq(den)));
                gradient[3] += main_term * (emax * (den - 1.0) * (ed50 - X[i]) / (Sq(delta * den)));
            }              

        }

        // Count endpoint
        if (endpoint_index == 3) {

            n_parameters = 4;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;   

            // Restricted estimates
            e0 = beta[0];    
            emax = beta[1];    
            ed50 = beta[2];    
            delta = beta[3];    

            for(i = 0; i < n; i++) {
                den = 1.0 + exp((ed50 - X[i]) / delta);
                mu = e0 + emax / den;
                if (mu <= 0.0) {
                    neg_log_likelihood += (lgamma(theta_vector[i]) + LogFactorial(Y[i]) - lgamma(theta_vector[i] + Y[i]) + (theta_vector[i] + Y[i]) * log(theta_vector[i] + exp(mu)) - Y[i] * mu - theta_vector[i] * log(theta_vector[i]));
                    main_term = theta_vector[i] * (exp(mu) - Y[i]) / (theta_vector[i] + exp(mu));
                } else {
                    neg_log_likelihood += (lgamma(theta_vector[i]) + LogFactorial(Y[i]) - lgamma(theta_vector[i] + Y[i]) + (theta_vector[i] + Y[i]) * (mu + log(theta_vector[i] * exp(-mu) + 1.0)) - Y[i] * mu - theta_vector[i] * log(theta_vector[i]));
                    main_term = theta_vector[i] * (1.0 - Y[i] * exp(-mu)) / (theta_vector[i] * exp(-mu) + 1.0);                    
                }
                gradient[0] += main_term;
                gradient[1] += main_term / den;
                gradient[2] += main_term * (- emax * (den - 1.0) / (delta * Sq(den)));
                gradient[3] += main_term * (emax * (den - 1.0) * (ed50 - X[i]) / (Sq(delta * den)));
            }
        }

        final_gradient = 0.0;
        for (i = 0; i < n_parameters; i++) final_gradient += abs(gradient[i]);

        const double f = neg_log_likelihood;
        grad.noalias() = gradient;

        return f;
    }
};

class RegressionSigEmax: public MFuncGrad
{
private:
    const MapVec X;
    const MapVec Y;
public:
    RegressionSigEmax(const MapVec x_, const MapVec y_) : X(x_), Y(y_) {}

    double f_grad(Constvec& beta, Refvec grad)
    {

        int i, n = Y.size(), n_parameters = 5;
        double neg_log_likelihood = 0.0, p, mu, e0, emax, ed50, h, sigma, den, main_term;
        Eigen::VectorXd gradient; 

        // Normal endpoint
        if (endpoint_index == 1) {

            n_parameters = 5;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;   

            // Restricted estimates
            e0 = beta[0];    
            emax = beta[1];    
            ed50 = max(beta[2], 0.01);    
            h = max(beta[3], 0.01);    
            sigma = max(beta[4], 0.0001);    

            for(i = 0; i < n; i++) {
                den = pow(X[i], h) + pow(ed50, h);
                mu = e0 + emax * pow(X[i], h) / den;
                neg_log_likelihood += (log(sqrt(2 * pi) * sigma) + Sq(Y[i] - mu) / (2 * Sq(sigma)));
                main_term = (mu - Y[i]) / Sq(sigma);
                gradient[0] += main_term;
                gradient[1] += main_term * pow(X[i], h) / den;
                gradient[2] += main_term * (- emax * h * pow(ed50, h - 1.0) * pow(X[i], h) / Sq(den));
                gradient[3] += main_term * (emax * pow(X[i], h) * pow(ed50, h) * (TruncLog(X[i] / ed50)) / Sq(den));
                gradient[4] += (1.0 / sigma - Sq(Y[i] - mu) / (sigma * sigma * sigma));
            }

        }

        // Binary endpoint
        if (endpoint_index == 2) {

            n_parameters = 4;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;           

            // Restricted estimates
            e0 = beta[0];    
            emax = beta[1];    
            ed50 = max(beta[2], 0.01);    
            h = max(beta[3], 0.01);    

            for(i = 0; i < n; i++) {
                den = pow(X[i], h) + pow(ed50, h);
                mu = e0 + emax * pow(X[i], h) / den;
                if (mu <= 0.0) {
                    p = exp(mu) / (1.0 + exp(mu));
                    neg_log_likelihood += (log(1.0 + exp(mu)) - Y[i] * mu);
                } else {
                    p = 1.0 / (1.0 + exp(-mu));
                    neg_log_likelihood += (mu + log(1.0 + exp(-mu)) - Y[i] * mu);
                } 
                main_term = p - Y[i];
                gradient[0] += main_term;
                gradient[1] += main_term * pow(X[i], h) / den;
                gradient[2] += main_term * (- emax * h * pow(ed50, h - 1.0) * pow(X[i], h) / Sq(den));
                gradient[3] += main_term * (emax * pow(X[i], h) * pow(ed50, h) * (TruncLog(X[i] / ed50)) / Sq(den));
            }              

        }

        // Count endpoint
        if (endpoint_index == 3) {

            n_parameters = 4;
            gradient.resize(n_parameters); 
            for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;   

            // Restricted estimates
            e0 = beta[0];    
            emax = beta[1];    
            ed50 = max(beta[2], 0.01);    
            h = max(beta[3], 0.01);    

            for(i = 0; i < n; i++) {
                den = pow(X[i], h) + pow(ed50, h);
                mu = e0 + emax * pow(X[i], h) / den;
                if (mu <= 0.0) {
                    neg_log_likelihood += (lgamma(theta_vector[i]) + LogFactorial(Y[i]) - lgamma(theta_vector[i] + Y[i]) + (theta_vector[i] + Y[i]) * log(theta_vector[i] + exp(mu)) - Y[i] * mu - theta_vector[i] * log(theta_vector[i]));
                    main_term = theta_vector[i] * (exp(mu) - Y[i]) / (theta_vector[i] + exp(mu));
                } else {
                    neg_log_likelihood += (lgamma(theta_vector[i]) + LogFactorial(Y[i]) - lgamma(theta_vector[i] + Y[i]) + (theta_vector[i] + Y[i]) * (mu + log(theta_vector[i] * exp(-mu) + 1.0)) - Y[i] * mu - theta_vector[i] * log(theta_vector[i]));
                    main_term = theta_vector[i] * (1.0 - Y[i] * exp(-mu)) / (theta_vector[i] * exp(-mu) + 1.0);                    
                }
                gradient[0] += main_term;
                gradient[1] += main_term * pow(X[i], h) / den;
                gradient[2] += main_term * (- emax * h * pow(ed50, h - 1.0) * pow(X[i], h) / Sq(den));
                gradient[3] += main_term * (emax * pow(X[i], h) * pow(ed50, h) * (TruncLog(X[i] / ed50)) / Sq(den));
            }

        }

        final_gradient = 0.0;
        for (i = 0; i < n_parameters; i++) final_gradient += abs(gradient[i]);

        const double f = neg_log_likelihood;
        grad.noalias() = gradient;

        return f;
    }
};

void SetInitialValues(vector<ModelInformation> &model_information, const vector<double> &dose, const vector<double> &response, const double &max_dose, const vector<int> &selected_models) {

    int i;
    vector<double> linear_coef(2);

    // Fit a linear regression model first    
    linear_coef = FitLinearModel(dose, response);

    double placebo_effect = linear_coef[0], max_effect = linear_coef[0] + linear_coef[1] * max_dose, d, d1, d2, p, p1, p2;

    // Initial values of the non-linear model parameters  
    vector<double> non_linear_parameters(2);

    ModelInformation current_model_information;

    for (i = 0; i < n_models; i++) {

        current_model_information.flag = selected_models[i];

        // Set initial values of the model parameters

        // Linear
        if (i == 0) {
            current_model_information.initial_values = linear_coef;
            current_model_information.n_parameters = 2;
        }

        // Quadratic
        if (i == 1) {
            non_linear_parameters[0] = 1.0;
            non_linear_parameters[1] = 0.0;
            current_model_information.initial_values = ComputeDoseResponseFunctionParameters(2, placebo_effect, max_effect, max_dose, non_linear_parameters);
            current_model_information.n_parameters = 3;
        }

        // Exponential
        if (i == 2) {
            non_linear_parameters[0] = max_dose;
            non_linear_parameters[1] = 0.0;
            current_model_information.initial_values = ComputeDoseResponseFunctionParameters(3, placebo_effect, max_effect, max_dose, non_linear_parameters);
            current_model_information.n_parameters = 3;
        }

        // Emax
        if (i == 3) {
            d = 0.5 * max_dose;
            p = 0.5;
            non_linear_parameters[0] = d * (1.0 - p) / p;
            non_linear_parameters[1] = 0.0;
            current_model_information.initial_values = ComputeDoseResponseFunctionParameters(4, placebo_effect, max_effect, max_dose, non_linear_parameters);
            current_model_information.n_parameters = 3;
        }

        // Logistic
        if (i == 4) {
            d1 = 0.33 * max_dose;
            p1 = 0.33;
            d2 = 0.66 * max_dose;
            p2 = 0.66;
            non_linear_parameters[0] = (d1 * Logit(p2) - d2 * Logit(p1)) / (Logit(p2) - Logit(p1));
            non_linear_parameters[1] = (d2 - d1)/(Logit(p2) - Logit(p1));
            current_model_information.initial_values = ComputeDoseResponseFunctionParameters(5, placebo_effect, max_effect, max_dose, non_linear_parameters);
            current_model_information.n_parameters = 4;
        }

        // SigEmax
        if (i == 5) {
            d1 = 0.33 * max_dose;
            p1 = 0.33;
            d2 = 0.66 * max_dose;
            p2 = 0.66;
            non_linear_parameters[1] = log(p1 * (1 - p2)/(p2 * (1 - p1)))/(log(d1/d2));
            non_linear_parameters[0] = pow((1 - p1)/p1, 1/non_linear_parameters[1]) * d1;
            current_model_information.initial_values = ComputeDoseResponseFunctionParameters(6, placebo_effect, max_effect, max_dose, non_linear_parameters);
            current_model_information.n_parameters = 4;
        }

        // Add a parameter that represents the standard deviation
        if (endpoint_index == 1) {
            current_model_information.n_parameters += 1;
            current_model_information.coef = FillVec(current_model_information.n_parameters, 0.0);
        } else {
            current_model_information.coef = FillVec(current_model_information.n_parameters, 0.0);
        }

        current_model_information.status = -1;
        current_model_information.criterion = 0.0;    
        current_model_information.target_dose = -1.0;
        current_model_information.convergence_criterion = -1.0;

        model_information[i] = current_model_information;        

    }

}


void FitDoseResponseModels(vector<ModelInformation> &model_information, const NumericVector &x, const NumericVector &y, const double &delta, const double &direction_index, const int &maxit) {

    const MapVec dose = Rcpp::as<MapVec>(x);
    const MapVec outcome = Rcpp::as<MapVec>(y);

    int i, j, n_parameters, convergence;

    /*******************************************************************************/

    // Convergence parameters
    // eps_f: Iteration stops if |f-f'|/|f|<eps_f, where f and f' are the current and previous value of the objective function (negative log likelihood) respectively.
    // eps_g: Iteration stops if ||g|| < eps_g * max(1, ||beta||), where beta is the current coefficient vector and g is the gradient.
    // maxit: Maximum number of iterations.

    double fopt, eps_f = 1e-08, eps_g = 1e-06;

    /*******************************************************************************/

    for (i = 0; i < n_models; i++) {

        convergence = 1;
        model_information[i].status = -1;
        model_information[i].convergence_criterion = -1.0;
        model_information[i].criterion = -1.0;
        model_information[i].target_dose = -1.0;

        if (model_information[i].flag == 1) {

            n_parameters = model_information[i].n_parameters;    
        
            Eigen::VectorXd beta(n_parameters);
            vector<double> coef(n_parameters);

            // Set the initial values
            if (endpoint_index == 1) {
                for (j = 0; j < n_parameters - 1; j++) beta[j] = model_information[i].initial_values[j];
                // Set the the initial value for the standard deviation
                beta[n_parameters - 1] = 1.0;    
            } else {
                for (j = 0; j < n_parameters; j++) beta[j] = model_information[i].initial_values[j];               
            }

            if (i == 0) {
                RegressionLinear linear_regression(dose, outcome);
                model_information[i].status = optim_lbfgs(linear_regression, beta, fopt, maxit, eps_f, eps_g);
            }

            if (i == 1) {
                iter = 1;
                RegressionQuadratic quadratic_regression(dose, outcome);
                model_information[i].status = optim_lbfgs(quadratic_regression, beta, fopt, maxit, eps_f, eps_g);
            }

            if (i == 2) {
                RegressionExponential exponential_regression(dose, outcome);
                model_information[i].status = optim_lbfgs(exponential_regression, beta, fopt, maxit, eps_f, eps_g);
                beta[2] = max(beta[2], 0.01);                                
            }

            if (i == 3) {
                RegressionEmax emax_regression(dose, outcome);
                model_information[i].status = optim_lbfgs(emax_regression, beta, fopt, maxit, eps_f, eps_g);
                beta[2] = max(beta[2], 0.01);                                
            }

            if (i == 4) {
                RegressionLogistic logistic_regression(dose, outcome);
                model_information[i].status = optim_lbfgs(logistic_regression, beta, fopt, maxit, eps_f, eps_g);
                beta[2] = max(beta[2], 0.01);             
                beta[3] = max(beta[3], 0.01);               
            }

            if (i == 5) {
                RegressionSigEmax sigemax_regression(dose, outcome);
                model_information[i].status = optim_lbfgs(sigemax_regression, beta, fopt, maxit, eps_f, eps_g);
                beta[2] = max(beta[2], 0.01);                
                beta[3] = max(beta[3], 0.01);               
            }

            if(model_information[i].status >= 0) {

                for (j = 0; j < n_parameters; j++) coef[j] = beta[j];

            }

            // Criteria for determining convergence
            if (isnan(fopt) || isnan(final_gradient) || abs(final_gradient) > 100.0 || model_information[i].status < 0) {
                convergence = 0;
                model_information[i].status = -1;
            }

            // Lack of convergence
            if (convergence == 1) {

                model_information[i].coef = coef;
                model_information[i].convergence_criterion = final_gradient;

                // Compute AIC
                model_information[i].criterion = 2.0 * fopt + 2.0 * (n_parameters + 0.0);
                model_information[i].target_dose = FindTargetDose(i + 1, coef, delta, direction_index);     
            }  

        }

    }
}

// [[Rcpp::export]]
List MCPModFitDRModels(const int &endpoint_index_arg, 
                       const IntegerVector &selected_models_arg, 
                       const NumericVector &x_arg, 
                       const NumericVector &y_arg, 
                       const double &delta_arg, 
                       const double &direction_index_arg,
                       const int &maxit_arg,
                       const NumericVector &theta_vector_arg) {


    // Global variable
    int endpoint_index_temp(endpoint_index_arg);
    endpoint_index = endpoint_index_temp;

    NumericVector x(x_arg);
    NumericVector y(y_arg);

    double delta(delta_arg);
    double direction_index(direction_index_arg);

    vector<int> selected_models = as<vector<int>>(selected_models_arg);

    int maxit(maxit_arg);

    theta_vector = as<vector<double>>(theta_vector_arg);

    /*******************************************************************************/

    // Compute initial values of the model parameters

    // Vector of model information parameters
    vector<ModelInformation> model_information(n_models);

    int i;
    vector<double> dose, response;
    dose = ConvertToDoubleVector(x);
    response = ConvertToDoubleVector(y);
    double max_dose = MaxVec(dose);

    SetInitialValues(model_information, dose, response, max_dose, selected_models);

    /*******************************************************************************/

    // Fit the dose-response models

    FitDoseResponseModels(model_information, x, y, delta, direction_index, maxit);

    /*******************************************************************************/

    // Save all results

    i = 0;
    Rcpp::List linear = List::create(                       
                       Named("status") = model_information[i].status,
                       Named("model_index") = i + 1,
                       Named("coef") = model_information[i].coef,
                       Named("criterion") = model_information[i].criterion,
                       Named("target_dose") = model_information[i].target_dose,
                       Named("initial_values") = model_information[i].initial_values,
                       Named("convergence_criterion") = model_information[i].convergence_criterion
                       );
    i = 1;
    Rcpp::List quadratic = List::create(                       
                       Named("status") = model_information[i].status,
                       Named("model_index") = i + 1,
                       Named("coef") = model_information[i].coef,
                       Named("criterion") = model_information[i].criterion,
                       Named("target_dose") = model_information[i].target_dose,
                       Named("initial_values") = model_information[i].initial_values,
                       Named("convergence_criterion") = model_information[i].convergence_criterion
                       );
    i = 2;
    Rcpp::List exponential = List::create(                       
                       Named("status") = model_information[i].status,
                       Named("model_index") = i + 1,
                       Named("coef") = model_information[i].coef,
                       Named("criterion") = model_information[i].criterion,
                       Named("target_dose") = model_information[i].target_dose,
                       Named("initial_values") = model_information[i].initial_values,
                       Named("convergence_criterion") = model_information[i].convergence_criterion
                       );
    i = 3;
    Rcpp::List emax = List::create(                       
                       Named("status") = model_information[i].status,
                       Named("model_index") = i + 1,
                       Named("coef") = model_information[i].coef,
                       Named("criterion") = model_information[i].criterion,
                       Named("target_dose") = model_information[i].target_dose,
                       Named("initial_values") = model_information[i].initial_values,
                       Named("convergence_criterion") = model_information[i].convergence_criterion
                       );
    i = 4;
    Rcpp::List logistic = List::create(                       
                       Named("status") = model_information[i].status,
                       Named("model_index") = i + 1,
                       Named("coef") = model_information[i].coef,
                       Named("criterion") = model_information[i].criterion,
                       Named("target_dose") = model_information[i].target_dose,
                       Named("initial_values") = model_information[i].initial_values,
                       Named("convergence_criterion") = model_information[i].convergence_criterion
                       );
    i = 5;
    Rcpp::List sigemax = List::create(                       
                       Named("status") = model_information[i].status,
                       Named("model_index") = i + 1,
                       Named("coef") = model_information[i].coef,
                       Named("criterion") = model_information[i].criterion,
                       Named("target_dose") = model_information[i].target_dose,
                       Named("initial_values") = model_information[i].initial_values,
                       Named("convergence_criterion") = model_information[i].convergence_criterion
                       );

    return List::create(
                        Named("linear") = linear,
                        Named("quadratic") = quadratic,
                        Named("exponential") = exponential,
                        Named("emax") = emax,
                        Named("logistic") = logistic,
                        Named("sigemax") = sigemax
                        );

    
}

// [[Rcpp::export]]
List MCPModRunSimulations(const int &endpoint_index_arg, 
                 const IntegerVector &selected_models_arg, 
                 const NumericVector &theta_arg, 
                 const NumericVector &theta_vector_arg, 
                 const double &delta_arg, 
                 const int &model_selection_index_arg, 
                 const List &opt_contrast_arg, 
                 const NumericVector &crit_value_arg, 
                 const List &sim_parameter_list_arg, 
                 const List &sim_model_list_arg, 
                 const double &direction_index_arg,
                 const double &go_threshold_arg,
                 const int &n_points_arg,
                 const int &maxit_arg) {

    // Global variable
    int endpoint_index_temp(endpoint_index_arg);
    endpoint_index = endpoint_index_temp;

    double delta(delta_arg);
    double direction_index(direction_index_arg);

    int model_selection_index(model_selection_index_arg);

    vector<int> selected_models = as<vector<int>>(selected_models_arg);

    Rcpp::List sim_parameter_list(sim_parameter_list_arg);
    Rcpp::List sim_model_list(sim_model_list_arg);
    Rcpp::List opt_contrast_list(opt_contrast_arg);

    // Vector of overdispersion parameters in the negative binomial distribution
    vector<double> theta = as<vector<double>>(theta_arg);

    // Long vector of overdispersion parameters in the negative binomial distribution (one value for each patient)
    theta_vector = as<vector<double>>(theta_vector_arg);

    // Multiplicity-adjusted critical value
    NumericVector crit_value(crit_value_arg);

    // Simulation parameters
    vector<int> n = as<vector<int>>(sim_parameter_list["n"]);
    vector<double> dose_levels = as<vector<double>>(sim_parameter_list["doses"]);
    vector<double> sd = as<vector<double>>(sim_model_list["sd"]);
    double dropout_rate = as<double>(sim_parameter_list["dropout_rate"]);
    int nsims = as<int>(sim_parameter_list["nsims"]);

    // Simulation models
    vector<int> sim_model_index = as<vector<int>>(sim_model_list["sim_model_index"]);
    NumericMatrix sim_parameter_values = sim_model_list["sim_parameter_values"];

    double go_threshold(go_threshold_arg);
    int n_points(n_points_arg);
    int n_scenarios = sim_parameter_values.nrow();
    int n_parameters = sim_parameter_values.ncol();

    double max_dose = dose_levels.back();

    int maxit(maxit_arg);

    /*******************************************************************************/

    // Other variables

    int i, j, k, index, sim, scenario, n_doses = dose_levels.size(), total_patients, n_selected_models = 0;
    vector<int> selected_index;

    for (i = 0; i < n_models; i++) {
        if (selected_models[i] == 1) {
            n_selected_models++;
            selected_index.push_back(i);
        }
    }

    vector<int> any_significant_models(nsims);
    vector<double> sign_model(n_selected_models), test_stat(n_selected_models), group_mean(n_doses), group_n(n_doses), var_group(n_doses), power(n_scenarios), beta(n_parameters), sample, model_weight(n_selected_models), temp_vec, true_target_dose(n_scenarios), go_prob(n_scenarios), dose, response;
    double denom, numer, pooled_variance, pred, temp_min, temp_max, mean, current_value, current_target_dose, group_logit, group_var, denominator, current_criterion, max_effect;

    NumericVector x, y;

    IntegerMatrix model_index(nsims, n_scenarios);

    NumericMatrix target_dose(nsims, n_scenarios), model_index_summary(n_scenarios, n_selected_models + 1), target_dose_summary(n_scenarios, 3), target_dose_categorical_summary(n_scenarios, n_doses + 1), dose_response_mean(n_points, n_scenarios), dose_response_lower(n_points, n_scenarios), dose_response_upper(n_points, n_scenarios);

    vector<ModelParameters> model_parameters;
    ModelParameters current_model_parameters;

    // X values for computing the estimated dose-response function
    vector<double> dosex(n_points); 

    for (i = 0; i < n_points; i++) {

        dosex[i] = max_dose * (i + 0.0) / (n_points - 1.0);     

    }

    // Set of optimal contrast matrices
    vector<OptimalContrastMatrix> optimal_contrast_matrix(n_scenarios);   

    // Populate optimal contrast matrices
    if (endpoint_index == 1) {
        NumericMatrix temp_matrix = opt_contrast_list[0];        
        optimal_contrast_matrix[0].matrix = temp_matrix;
    } else {
        for (i = 0; i < n_scenarios; i++) {
            NumericMatrix temp_matrix = opt_contrast_list[i];        
            optimal_contrast_matrix[i].matrix = temp_matrix;
        }
    }

    /*******************************************************************************/

    // Vector of model information parameters

    vector<ModelInformation> model_information(n_models);

    /*******************************************************************************/

    for (scenario = 0; scenario < n_scenarios; scenario++) {

        // Scenario-specific model parameters
        for (i = 0; i < n_parameters; i++) beta[i] = sim_parameter_values(scenario, i);

        // Find the true target dose    
        true_target_dose[scenario] = FindTargetDose(sim_model_index[0], beta, delta, direction_index);

        model_parameters.clear();

        // Run simulations

        for (sim = 0; sim < nsims; sim++) {

            // Generate a data set

            x.erase(x.begin(), x.end());
            y.erase(y.begin(), y.end());

            total_patients = 0;
            current_value = 0;

            for (j = 0; j < n_doses; j++) {

                group_mean[j] = 0.0;
                group_n[j] = 0.0;
                mean = DoseResponseFunction(dose_levels[j], sim_model_index[0], beta, direction_index);

                for (k = 0; k < n[j]; k++) {

                    // Remove dropouts
                    if (rcpp_uniform(1, 0.0, 1.0)[0] >= dropout_rate) {
                        x.push_back(dose_levels[j]);

                        // Normal endpoint
                        if (endpoint_index == 1) {
                            current_value = rcpp_normal(1, mean, sd[j])[0];
                            y.push_back(current_value);
                        }

                        // Binary endpoint
                        if (endpoint_index == 2) {
                            current_value = rcpp_binary(1, mean)[0];
                            y.push_back(current_value);
                        }

                        // Count endpoint
                        if (endpoint_index == 3) {
                            current_value = rcpp_nbinom(1, theta[j], theta[j] /(theta[j] + mean))[0]; 
                            y.push_back(current_value);
                        }

                        group_mean[j] += current_value;
                        group_n[j]++;
                        total_patients++;
                    }

                }

                // Compute the group-specific mean
                if (group_n[j] > 0.0) group_mean[j] /= group_n[j]; 

            }


            // Normal endpoint
            pooled_variance = 0.0;
            if (endpoint_index == 1) {
                for (j = 0; j < n_doses; j++) {
                  for (i = 0; i < total_patients; i++) {
                    if (x[i] == dose_levels[j]) {
                      pooled_variance += Sq(y[i] - group_mean[j]);     
                    }
                  }
                }
                if ((double) (total_patients - n_doses) > 0.0) pooled_variance /= ((double) (total_patients - n_doses));
            }

            // Compute optimal contrasts
            group_logit = 0.0;
            group_var = 0.0;

            // Normal endpoint
            if (endpoint_index == 1) {

                for (i = 0; i < n_selected_models; i++) {
                  denom = 0.0;
                  numer = 0.0;
                  for (j = 0; j < n_doses; j++) {
                      numer = numer + optimal_contrast_matrix[0].matrix(j, i) * group_mean[j];
                      if (group_n[j] > 0.0) denom = denom + Sq(optimal_contrast_matrix[0].matrix(j, i)) / group_n[j];
                  }
                  test_stat[i] = numer / sqrt(pooled_variance * denom);
                  if (direction_index == 1) sign_model[i] = (test_stat[i] >= crit_value[0]);
                  if (direction_index == -1) sign_model[i] = (test_stat[i] <= crit_value[0]);
                }

            }


            // Binary endpoint
            if (endpoint_index == 2) {

                for (i = 0; i < n_selected_models; i++) {
                  denom = 0.0;
                  numer = 0.0;
                  for (j = 0; j < n_doses; j++) {
                      if (group_mean[j] == 0) group_mean[j] = 1.0 / (3.0 * group_n[j] + 2.0);
                      if (group_mean[j] == group_n[j]) group_mean[j] = (3.0 * group_n[j] + 1.0) / (3.0 * group_n[j] + 2.0);
                      if (group_mean[j] / (1.0 - group_mean[j]) > 0.0) group_logit = log(group_mean[j] / (1.0 - group_mean[j]));
                      if (group_n[j] > 0.0 && group_mean[j] > 0.0) group_var = 1.0 / (group_n[j] * group_mean[j] * (1 - group_mean[j]));
                      numer = numer + optimal_contrast_matrix[scenario].matrix(j, i) * group_logit;
                      denom = denom + Sq(optimal_contrast_matrix[scenario].matrix(j, i)) * group_var;
                  }
                  test_stat[i] = numer / sqrt(denom);
                  if (direction_index == 1) sign_model[i] = (test_stat[i] >= crit_value[scenario]);
                  if (direction_index == -1) sign_model[i] = (test_stat[i] <= crit_value[scenario]);
                }

            }

            // Count endpoint
            if (endpoint_index == 3) {

                for (i = 0; i < n_selected_models; i++) {
                  denom = 0.0;
                  numer = 0.0;
                  for (j = 0; j < n_doses; j++) {
                      if (group_mean[j] == 0) group_mean[j] = rcpp_qgamma(0.5, 1.0 / 3.0, 1.0 / group_n[j]);
                      if (group_n[j] > 0.0 && group_mean[j] > 0.0) group_var = (theta[j] + group_mean[j]) / (group_n[j] * theta[j] * group_mean[j]);
                      numer = numer + optimal_contrast_matrix[scenario].matrix(j, i) * log(group_mean[j]);
                      denom = denom + Sq(optimal_contrast_matrix[scenario].matrix(j, i)) * group_var;
                  }
                  test_stat[i] = numer / sqrt(denom);
                  if (direction_index == 1) sign_model[i] = (test_stat[i] >= crit_value[scenario]);
                  if (direction_index == -1) sign_model[i] = (test_stat[i] <= crit_value[scenario]);
                }

            }

            // Find if there are any significant models
            any_significant_models[sim]= 0;
            for (i = 0; i < n_selected_models; i++) {
                if (sign_model[i] == 1) any_significant_models[sim] = 1;
            }

            /*******************************************************************************/

            // Fit the dose-response models 

            // Model selection and target dose estimation
            model_index(sim, scenario) = -1;
            target_dose(sim, scenario) = -1.0;

            if (any_significant_models[sim] == 1) {

                power[scenario]++;

                // Compute initial values of the model parameters

                // Vector of model information parameters
                vector<ModelInformation> model_information(n_models);

                dose = ConvertToDoubleVector(x);
                response = ConvertToDoubleVector(y);

                SetInitialValues(model_information, dose, response, max_dose, selected_models);

                FitDoseResponseModels(model_information, x, y, delta, direction_index, maxit);

                // Find the significant model with the smallest criterion
                if (model_selection_index == 1) {

                    temp_min = 0.0;
    
                    for (i = 0; i < n_selected_models; i++) {

                        if (model_index(sim, scenario) == -1 && sign_model[i] == 1 && model_information[selected_index[i]].status >= 0) {
                            temp_min = model_information[selected_index[i]].criterion;
                            model_index(sim, scenario) = i + 1;
                            target_dose(sim, scenario) = model_information[selected_index[i]].target_dose;
                        }

                        if (model_index(sim, scenario) >= 0 && sign_model[i] == 1 && model_information[selected_index[i]].status >= 0) {
                            if (temp_min > model_information[selected_index[i]].criterion) {
                                model_index(sim, scenario) = i + 1;
                                temp_min = model_information[selected_index[i]].criterion;
                                target_dose(sim, scenario) = model_information[selected_index[i]].target_dose;
                            }
                        }
                    }

                }

                // Find the significant model with the largest test statistic
                if (model_selection_index == 2) {

                    temp_max = 0.0;
    
                    for (i = 0; i < n_selected_models; i++) {

                        if (model_index(sim, scenario) == -1 && sign_model[i] == 1 && model_information[selected_index[i]].status >= 0) {
                            temp_max = test_stat[i];
                            model_index(sim, scenario) = i + 1;
                            target_dose(sim, scenario) = model_information[selected_index[i]].target_dose;
                        }

                        if (model_index(sim, scenario) >= 0 && sign_model[i] == 1 && model_information[selected_index[i]].status >= 0) {
                            if (temp_max < test_stat[i]) {
                                model_index(sim, scenario) = i + 1;
                                temp_max = test_stat[i];
                                target_dose(sim, scenario) = model_information[selected_index[i]].target_dose;
                            }
                        }
                    }

                }


                // Find the target dose based on the weighted criterion
                if (model_selection_index == 3) {

                    target_dose(sim, scenario) = 0.0;

                    for (i = 0; i < n_selected_models; i++) {

                        model_weight[i] = 0.0;

                        if (sign_model[i] == 1 && model_information[selected_index[i]].status >= 0) {

                            current_criterion = model_information[selected_index[i]].criterion;
                            denominator = 0.0;

                            for (j = 0; j < n_selected_models; j++) {

                                if (sign_model[j] == 1 && model_information[selected_index[j]].status >= 0) denominator += exp(- 0.5 * (model_information[selected_index[j]].criterion - current_criterion));

                            }

                            if (abs(denominator) > 0.0001) model_weight[i] = 1.0 / denominator;

                            target_dose(sim, scenario) += (model_information[selected_index[i]].target_dose * model_weight[i]);

                        }
                   

                    }

                    if (target_dose(sim, scenario) == 0) target_dose(sim, scenario) = -1.0;

                    model_index(sim, scenario) = 1;

                }

                // Save the model parameters of the best dose-response model
                if (model_selection_index != 3) {

                    index = model_index(sim, scenario);

                    if (index >= 1) {
                        current_model_parameters.scenario = scenario;
                        current_model_parameters.model = selected_index[index - 1] + 1;
                        current_model_parameters.coef = model_information[selected_index[index - 1]].coef;
                        model_parameters.push_back(current_model_parameters);

                        // Maximum effect for the best dose-response curve
                        max_effect = FindMaxEffect(current_model_parameters.model, current_model_parameters.coef, max_dose, direction_index);

                        if (direction_index == 1) {
                            if (max_effect >= go_threshold) go_prob[scenario]++; 
                        }    

                        if (direction_index == -1) {
                            if (max_effect <= go_threshold) go_prob[scenario]++; 
                        }                     

                    }

                }

            }

            // Summary of model selection
            if (any_significant_models[sim] == 1 && model_selection_index != 3 && model_index(sim, scenario) >= 1) {

                model_index_summary(scenario, model_index(sim, scenario))++;                               

            } else {

                model_index_summary(scenario, 0)++;                               

            } 

        }
        // End of the simulation loop

        for (i = 0; i < 3; i++) {
            target_dose_summary(scenario, i) = NA_REAL;
        }

        for (i = 0; i < n_points; i++) {
            dose_response_lower(i, scenario) = NA_REAL;
            dose_response_mean(i, scenario) = NA_REAL;
            dose_response_upper(i, scenario) = NA_REAL;            
        }

        for (i = 0; i < n_doses + 1; i++) target_dose_categorical_summary(scenario, i) = 0.0;

        // Summary of target dose estimates
        temp_vec.clear();

        // Extract non-missing value of the target dose
        for (sim = 0; sim < nsims; sim++) {

            current_target_dose = target_dose(sim, scenario);

            if (any_significant_models[sim] == 1 && current_target_dose >= 0.0) {

                // Data for a continuous summary of target doses
                temp_vec.push_back(current_target_dose);

                // Data for a categorical summary of target doses

                // Target dose is greater than the largest dose
                if (current_target_dose > max_dose) target_dose_categorical_summary(scenario, n_doses)++;

                // Target dose is within the dose range
                if (current_target_dose <= max_dose) {

                        // Find the dose interval
                        for (j = 0; j < n_doses - 1; j++) {

                            if (current_target_dose > dose_levels[j] && current_target_dose <= dose_levels[j + 1]) target_dose_categorical_summary(scenario, j + 1)++;

                        }

                }
            } else {

                // Target dose was not found
                target_dose_categorical_summary(scenario, 0)++; 

            }

        }

        if (temp_vec.size() > 2) {

            // Lower 95% bound for the target dose
            target_dose_summary(scenario, 0) = Quantile(temp_vec, 0.025);

            // Mean target dose
            target_dose_summary(scenario, 1) = MeanVec(temp_vec);

            // Upper 95% bound for the target dose
            target_dose_summary(scenario, 2) = Quantile(temp_vec, 0.975);

        }

        // Compute the estimated dose-response function and confidence band for the current scenario
        if (model_parameters.size() > 0 && model_selection_index != 3) {

             for (i = 0; i < n_points; i++) {

                temp_vec.clear();

                for (j = 0; j < model_parameters.size(); j++) {

                    pred = DoseResponseFunction(dosex[i], model_parameters[j].model, model_parameters[j].coef, direction_index);

                    temp_vec.push_back(pred); 
 
                }

            if (temp_vec.size() > 0) {

                    // Compute the mean and 95% limits    
                    dose_response_lower(i, scenario) = Quantile(temp_vec, 0.025);
                    dose_response_mean(i, scenario) = MeanVec(temp_vec);
                    dose_response_upper(i, scenario) = Quantile(temp_vec, 0.975);
            } 



          }


        }
 
    }
    // End of the scenario loop

    // Compute summaries of simulation results
    power = AverageSimResultsVector(power, nsims);

    model_index_summary = AverageSimResultsMatrix(model_index_summary, nsims);

    target_dose_categorical_summary = AverageSimResultsMatrix(target_dose_categorical_summary, nsims);
    go_prob = AverageSimResultsVector(go_prob, nsims);

    return List::create(Named("power") = power,
                        Named("dosex") = dosex,
                        Named("go_prob") = go_prob,
                        Named("dose_response_upper") = dose_response_upper,
                        Named("dose_response_mean") = dose_response_mean,
                        Named("dose_response_lower") = dose_response_lower,
                        Named("model_index_summary") = model_index_summary,
                        Named("target_dose_summary") = target_dose_summary,
                        Named("target_dose_categorical_summary") = target_dose_categorical_summary,
                        Named("true_target_dose") = true_target_dose
                        );




}
// End of Simulation

/*********************************************/

// Test functions

// [[Rcpp::export]]
double TestDoseResponseFunction(const double &x, const double &model, const NumericVector &coef) {

    int i, m = coef.size();
    vector<double> local_coef(m);

    for (i = 0; i < m; i++) local_coef[i] = coef[i];

    endpoint_index = 1;

    double y = DoseResponseFunction(x, model, local_coef, 1); 

    return y;

}

// [[Rcpp::export]]
double TestFindTargetDose(const double &delta, const int &model, const NumericVector &coef) {

    int i, m = coef.size();
    vector<double> local_coef(m);

    for (i = 0; i < m; i++) local_coef[i] = coef[i];

    endpoint_index = 1;

    double dose = FindTargetDose(model, local_coef, delta, 1); 

    return dose;

}
