//
// Felix Rouleau; 20188552; felix.rouleau@gmail.com
//
#include <math.h>
#include <stdio.h>

enum derivativeTypes {
    QUESTION_A,
    QUESTION_B,
    QUESTION_C
};

/// \brief Approximate the cmv using Newton's method
/// \param cmv Current value of cmv
/// \param y Values of y
/// \param lenY Length of y
/// \return Value of the function at cmv
double newtonApproximation(double cmv, double *y, int lenY) {
    // Function is sum from 1->N (yi^cmv*ln(yi)) / sum 1->N yi^cmv - 1/cmv * "constant"
    // "Constant" is 1/N sum 1->N ln(yi).
    double total_num = 0;
    double total_denum = 0;
    double rightHand = 0;
    double total;
    for (int i = 0; i < lenY; i++) {
        total_num += pow(y[i], cmv) * log(y[i]);
        total_denum += pow(y[i], cmv);
        rightHand += log(y[i]);
    }
    total = total_num / total_denum;
    total -= 1 / cmv;
    total -= rightHand / lenY;
    return total;
}

/// \brief Derivative of the given function. Returns a different derivative based on the type parameter
/// \param cmv Current value of cmv
/// \param y Values of y
/// \param lenY Length of y
/// \param type Derivative type to use
/// \return Value of the derivative at cmv
double newtonDerivative(double cmv, double *y, int lenY, enum derivativeTypes type) {
    switch (type) {
        case QUESTION_A:
            // Version where it is f'(x) = [f(x+h) - f(x)] / h where h is 10e-6
            return (newtonApproximation(cmv + 10e-6, y, lenY) - newtonApproximation(cmv, y, lenY)) / 10e-6;
        case QUESTION_B:
            // Version where it is f'(x) = [-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)] / 12h where h is 10e-6
            return (-newtonApproximation(cmv + 2 * 10e-6, y, lenY) +
                    8 * newtonApproximation(cmv + 10e-6, y, lenY) -
                    8 * newtonApproximation(cmv - 10e-6, y, lenY) +
                    newtonApproximation(cmv - 2 * 10e-6, y, lenY)) / (12 * 10e-6);
        case QUESTION_C:
        {
            // Version where it is the analytic derivative of f(x)

            // First and Second Summation are left hand side
            // First Summation
            double total;
            double temp = 0;
            for (int i = 0; i < lenY; i++) {
                temp += pow(y[i], cmv) * pow(log(y[i]), 2);
            }
            total = temp;

            // Second Summation
            temp = 0;
            for (int i = 0; i < lenY; i++) {
                temp += pow(y[i], cmv);
            }
            total *= temp;

            // Third and Fourth Summation are right hand side
            // Third Summation
            temp = 0;
            double temp_right = 0;
            for (int i = 0; i < lenY; i++) {
                temp += pow(y[i], cmv) * log(y[i]);
            }
            temp_right = temp;

            // Fourth Summation
            temp = 0;
            for (int i = 0; i < lenY; i++) {
                temp += pow(y[i], cmv) * log(y[i]);
            }
            temp_right *= temp;
            total -= temp_right;

            // Fifth Summation (Bottom)
            temp = 0;
            for (int i = 0; i < lenY; i++) {
                temp += pow(y[i], cmv);
            }
            temp = pow(temp, 2);
            total /= temp;

            // Add 1/cmv^2
            total += 1 / pow(cmv, 2);

            return total;
        }

    }
}

/// \brief Estimate the cmv using Newton's method. Choose the derivative based on the type parameter
/// \param y Values of y
/// \param lenY Length of y
/// \param type Derivative type to use
/// \return Estimated cmv
double cmvEstimator(double *y, int lenY, enum derivativeTypes type) {
    double tolerance = 10e-6;
    double cmv = 0.25;
    double old_cmv;
    do {
        old_cmv = cmv;
        cmv = cmv - newtonApproximation(cmv, y, lenY) / newtonDerivative(cmv, y, lenY, type);
    } while (fabs(cmv - old_cmv) > tolerance &&
             newtonApproximation(cmv, y, lenY) > tolerance &&
             newtonDerivative(cmv, y, lenY, type) != 0);
}

int main() {
    double y[] = {0.11, 0.24, 0.27, 0.52, 1.13, 1.54, 1.71, 1.84, 1.92, 2.01};
    int lenY = 10;
    double cmv = cmvEstimator(y, lenY, QUESTION_A);
    printf("CMV for Question A is %f\n", cmv);
    cmv = cmvEstimator(y, lenY, QUESTION_B);
    printf("CMV for Question B is %f\n", cmv);
    cmv = cmvEstimator(y, lenY, QUESTION_C);
    printf("CMV for Question C is %f", cmv);
    return 0;
}
