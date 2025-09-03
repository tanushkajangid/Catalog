#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <iomanip>

using namespace std;

// Function to convert from any base to decimal
long long baseToDecimal(const string& value, int base) {
    long long result = 0;
    long long power = 1;
    
    for (int i = value.length() - 1; i >= 0; i--) {
        int digit;
        char c = value[i];
        
        if (c >= '0' && c <= '9') {
            digit = c - '0';
        } else if (c >= 'A' && c <= 'Z') {
            digit = c - 'A' + 10;
        } else if (c >= 'a' && c <= 'z') {
            digit = c - 'a' + 10;
        } else {
            continue;
        }
        
        if (digit >= base) continue;
        
        result += digit * power;
        power *= base;
    }
    return result;
}

// Gaussian elimination to solve system of equations
vector<double> gaussianElimination(vector<vector<double>>& matrix) {
    int n = matrix.size();
    
    // Forward elimination
    for (int i = 0; i < n; i++) {
        // Find pivot
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(matrix[k][i]) > abs(matrix[maxRow][i])) {
                maxRow = k;
            }
        }
        swap(matrix[maxRow], matrix[i]);
        
        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < n; k++) {
            double factor = matrix[k][i] / matrix[i][i];
            for (int j = i; j <= n; j++) {
                matrix[k][j] -= factor * matrix[i][j];
            }
        }
    }
    
    // Back substitution
    vector<double> solution(n);
    for (int i = n - 1; i >= 0; i--) {
        solution[i] = matrix[i][n];
        for (int j = i + 1; j < n; j++) {
            solution[i] -= matrix[i][j] * solution[j];
        }
        solution[i] /= matrix[i][i];
    }
    
    return solution;
}

// Method 1: Vandermonde Matrix Method
double solveVandermonde(vector<pair<int, long long>>& points) {
    int n = points.size();
    vector<vector<double>> matrix(n, vector<double>(n + 1));
    
    cout << "\n=== Vandermonde Matrix Method ===" << endl;
    cout << "Building system: [1 x x^2 ... x^(n-1)] * [a0 a1 a2 ... a(n-1)] = [y]" << endl;
    
    for (int i = 0; i < n; i++) {
        int x = points[i].first;
        long long y = points[i].second;
        
        // Fill row: [1, x, x^2, x^3, ..., x^(n-1), y]
        for (int j = 0; j < n; j++) {
            matrix[i][j] = pow(x, j);
        }
        matrix[i][n] = y; // RHS
        
        cout << "Row " << i << ": ";
        for (int j = 0; j < n; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << "= " << y << endl;
    }
    
    vector<double> coefficients = gaussianElimination(matrix);
    
    cout << "\nPolynomial coefficients:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "a" << i << " = " << coefficients[i] << endl;
    }
    
    return coefficients[0]; // Constant term is a0
}

// Method 2: Newton's Divided Differences
double newtonDividedDifference(vector<pair<int, long long>>& points) {
    int n = points.size();
    vector<vector<double>> table(n, vector<double>(n, 0));
    
    cout << "\n=== Newton's Divided Differences Method ===" << endl;
    
    // Fill first column with y values
    for (int i = 0; i < n; i++) {
        table[i][0] = points[i].second;
    }
    
    // Fill the divided differences table
    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / 
                         (points[i + j].first - points[i].first);
        }
    }
    
    cout << "Divided Differences Table:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n - i; j++) {
            cout << setw(10) << fixed << setprecision(2) << table[i][j] << " ";
        }
        cout << endl;
    }
    
    // Calculate constant term by evaluating at x=0
    double result = table[0][0]; // f[x0]
    double product = 1.0;
    
    for (int i = 1; i < n; i++) {
        product *= (0 - points[i - 1].first); // (0 - x0)(0 - x1)...(0 - x_{i-1})
        result += table[0][i] * product;
    }
    
    cout << "\nNewton form evaluation at x=0:" << endl;
    cout << "f(0) = " << table[0][0];
    product = 1.0;
    for (int i = 1; i < n; i++) {
        product *= (0 - points[i - 1].first);
        cout << " + " << table[0][i] << "*" << product;
    }
    cout << " = " << result << endl;
    
    return result;
}

// Method 3: Direct substitution verification
void verifyPolynomial(vector<pair<int, long long>>& points, vector<double>& coeffs) {
    cout << "\n=== Polynomial Verification ===" << endl;
    cout << "Polynomial: f(x) = ";
    for (int i = coeffs.size() - 1; i >= 0; i--) {
        if (i == coeffs.size() - 1) {
            cout << coeffs[i];
        } else {
            cout << " + " << coeffs[i];
        }
        if (i > 0) cout << "*x";
        if (i > 1) cout << "^" << i;
    }
    cout << endl;
    
    cout << "\nVerifying points:" << endl;
    for (auto& point : points) {
        int x = point.first;
        long long expected = point.second;
        
        double calculated = 0;
        for (int i = 0; i < coeffs.size(); i++) {
            calculated += coeffs[i] * pow(x, i);
        }
        
        cout << "f(" << x << ") = " << calculated 
             << ", expected = " << expected 
             << ", match = " << (abs(calculated - expected) < 1e-9 ? "YES" : "NO") << endl;
    }
}

int main() {
    cout << "Testing Multiple Methods for Test Case 2" << endl;
    cout << "=========================================" << endl;
    
    // Test Case 2 data
    map<int, pair<int, string>> roots = {
        {1, {10, "4"}},
        {2, {2, "111"}},
        {3, {10, "12"}},
        {6, {4, "213"}}
    };
    
    int k = 3; // Use first 3 points
    
    vector<pair<int, long long>> points;
    
    cout << "\nDecoding roots:" << endl;
    for (auto& root : roots) {
        int x = root.first;
        int base = root.second.first;
        string value = root.second.second;
        
        long long y = baseToDecimal(value, base);
        points.push_back({x, y});
        
        cout << "x=" << x << ", base=" << base << ", value=\"" << value 
             << "\" -> y=" << y << endl;
    }
    
    // Use first k points
    vector<pair<int, long long>> selectedPoints(points.begin(), points.begin() + k);
    
    cout << "\nUsing points: ";
    for (auto& p : selectedPoints) {
        cout << "(" << p.first << "," << p.second << ") ";
    }
    cout << endl;
    
    // Method 1: Vandermonde Matrix
    double constant1 = solveVandermonde(selectedPoints);
    
    // Method 2: Newton's Divided Differences  
    double constant2 = newtonDividedDifference(selectedPoints);
    
    cout << "\n" << string(50, '=') << endl;
    cout << "RESULTS COMPARISON:" << endl;
    cout << string(50, '=') << endl;
    cout << "Vandermonde Matrix Method: " << constant1 << endl;
    cout << "Newton's Method: " << constant2 << endl;
    cout << "Rounded constant: " << (long long)round(constant1) << endl;
    
    return 0;
}