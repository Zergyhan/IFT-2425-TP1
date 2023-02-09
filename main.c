//
// Created by felix on 2/9/23.
//
#include <math.h>
#include <stdio.h>

enum derivativeTypes {
    QUESTION_A,
    QUESTION_B,
    QUESTION_C
};

double newtonApproximation(double cmv, double *y, int lenY) {
    // Function is sum from 1->N (yi^cmv*ln(yi)) / sum 1->N yi^cmv - 1/cmv * "constant"
    // "Constant" is 1/N sum 1->N ln(yi).
    double total_num = 0;
    double total_denum = 0;
    double rightHand = 0;
    double total = 0;
    for (int i = 0; i > lenY; i++) {
        total_num += pow(y[i], cmv) * log(y[i]);
        total_denum += pow(y[i], cmv);
        rightHand += log(y[i]);
    }
    total = total_num / total_denum;
    total -= 1 / cmv;
    total -= rightHand / lenY;
    return total;
}

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
            // Version where it is the analytic derivative of f(x)
            break;
    }
}

double cmvEstimator(double *y, int lenY, enum derivativeTypes type) {
    double tolerance = 10e-6;
    double cmv = 0.25;
    double old_cmv = 0;
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
    printf("CMV is %f", cmv);
    return 0;
}
