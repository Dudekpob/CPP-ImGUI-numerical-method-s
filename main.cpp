// Dear ImGui: standalone example application for DirectX 9
// If you are new to Dear ImGui, read documentation from the docs/ folder + read the top of imgui.cpp.
// Read online: https://github.com/ocornut/imgui/tree/master/docs

#include "imgui.h"
#include "imgui_impl_dx9.h"
#include "imgui_impl_win32.h"
#include <d3d9.h>
#include <tchar.h>
#include <iostream>
#include "Stategy.h"
#include <iomanip>
#include "eigen3/Eigen/Core"
#include <vector>
#include "eigen3/Eigen/Cholesky"
#include "eigen3/Eigen/QR"
#include "eigen3/Eigen/SVD"
#include <random>
#include <math.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <string>
#include <cctype>
#include <iterator>
#include "tinyexpr.h"


#define M_PI 3.14159265358979323846


// Data
static LPDIRECT3D9              g_pD3D = NULL;
static LPDIRECT3DDEVICE9        g_pd3dDevice = NULL;
static D3DPRESENT_PARAMETERS    g_d3dpp = {};

const double goldenRatio = 0.618;
const double accuracy = 0.00001;
const double lowerBound = -1;
const double upperBound = 1;


// Forward declarations of helper functions
bool CreateDeviceD3D(HWND hWnd);
void CleanupDeviceD3D();
void ResetDevice();

LRESULT WINAPI WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);
char name[30] = "";
void fn() {
    ImGui::InputText("Equation", name, IM_ARRAYSIZE(name));
}

double solve = 0;

double root_bs(double(*f)(double), double a, double b, double eps, int& flag)
{
    double xl, x0, xr;
    int i, iter = 1000;           /* Max number of iterations */

/* check the bisection condition */
    if (f(a) * f(b) > 0.0)
    {
        flag = 0;
        return 0.0;
    }

    /* finding a single root */
    i = 0;
    xl = a;
    xr = b;

    while (fabs(xr - xl) >= eps)
    {
        i = i + 1;
        x0 = (xr + xl) / 2.0;
        if ((f(xl) * f(x0)) <= 0.0)
            xr = x0;
        else
            xl = x0;
        if (i >= iter) break;
    }
    flag = i;
    return x0;
}



double fzero(double x)
{
    std::cout << x << std::endl;
   


    const char* expression = name;
    double y;
    te_parser tep;
    tep.set_vars({ {"x", &x} });
    /* Compile the expression with variables. */
    if (tep.compile(expression)) {
         double y = tep.evaluate();
        int k = 0;
        return y;
    }
    
  
        
    
    return 0;
    
    //y = exp(x) * log(x) - cos(x * x);

    //    y = log(x*x+2.0)*cos(x) + sin(x);
    //    y = 2.0*x*x*x + 6.0*x -21;
    //    y = x*x - 6.0*x +9.0;
  
}


int bisection_caclulation()
{
    double a, b, root, eps;
    int flag;
    std::cout.setf(std::ios::fixed | std::ios::showpoint);
    std::cout.precision(5);

    a = 1.0;                    // left endpoint
    b = 4.0;                    // right endpoint
    eps = 1.0e-6;               // desired uncertainity of the root

    root = root_bs(fzero, a, b, eps, flag);
    solve = root;
    if (flag == 0) std::cout << " no root for Bisectional method" << std::endl;
    else
    {
        std::cout << " iterations" << "      root" << std::endl;
        std::cout << std::setw(7) << flag << std::setw(14) << root << std::endl;
    }
    return 0;
}



class CostFunction {
public:
    CostFunction(double* a, double* b, double* c, int max_iter, double min_step, bool is_out) :
        a_(a), b_(b), c_(c), max_iter_(max_iter), min_step_(min_step), is_out_(is_out)
    {}

    void addObservation(double x, double y)
    {
        std::vector<double> ob;
        ob.push_back(x);
        ob.push_back(y);
        obs_.push_back(ob);
    }

    void calcJ_fx()
    {
        J_.resize(obs_.size(), 3);
        fx_.resize(obs_.size(), 1);

        for (size_t i = 0; i < obs_.size(); i++)
        {
            std::vector<double>& ob = obs_.at(i);
            double& x = ob.at(0);
            double& y = ob.at(1);
            double j1 = -x * x * exp(*a_ * x * x + *b_ * x + *c_);
            double j2 = -x * exp(*a_ * x * x + *b_ * x + *c_);
            double j3 = -exp(*a_ * x * x + *b_ * x + *c_);
            J_(i, 0) = j1;
            J_(i, 1) = j2;
            J_(i, 2) = j3;
            fx_(i, 0) = y - exp(*a_ * x * x + *b_ * x + *c_);
        }
    }

    void calcH_b()
    {
        H_ = J_.transpose() * J_;
        B_ = -J_.transpose() * fx_;
    }

    void calcDeltax()
    {
        deltax_ = H_.ldlt().solve(B_);
    }

    void updateX()
    {
        *a_ += deltax_(0);
        *b_ += deltax_(1);
        *c_ += deltax_(2);
    }

    double getCost()
    {
        Eigen::MatrixXd cost = fx_.transpose() * fx_;
        return cost(0, 0);
    }

    void solveByGaussNewton()
    {
        double sumt = 0;
        bool is_conv = false;
        for (size_t i = 0; i < max_iter_; i++)
        {
            calcJ_fx();
            calcH_b();
            calcDeltax();
            double delta = deltax_.transpose() * deltax_;
            if (is_out_)
            {
                std::cout << "Iter: " << std::left << std::setw(3) << i << " Result: " << std::left << std::setw(10) << *a_ << " " << std::left << std::setw(10) << *b_ << " " << std::left << std::setw(10) << *c_ <<
                    " step: " << std::left << std::setw(14) << delta << " cost: " << std::left << std::setw(14) << getCost() << " time: " << std::left << std::setw(14) <<
                    " total_time: " << std::left << std::setw(14) << std::endl;
            }
            if (delta < min_step_)
            {
                is_conv = true;
                break;
            }
            updateX();
               solve = delta;
        }
     
        if (is_conv == true)
            std::cout << "\nConverged\n";
        else
            std::cout << "\nDiverged\n\n";
    }

    Eigen::MatrixXd fx_;
    Eigen::MatrixXd J_; // 雅克比矩阵
    Eigen::Matrix3d H_; // H矩阵
    Eigen::Vector3d B_;
    Eigen::Vector3d deltax_;
    std::vector< std::vector<double>  > obs_; // 观测
    double* a_, * b_, * c_;

    int max_iter_;
    double min_step_;
    bool is_out_;
};//class CostFunction



int newton() {

    const double aa = 0.1, bb = 0.5, cc = 2; 
    double a = 0.0, b = 0.0, c = 0.0; 

   
    CostFunction cost_func(&a, &b, &c, 50, 1e-10, true);
    const char* expression = name;
    double y;
    te_parser tep;
 
    const size_t N = 100; 
    for (size_t i = 0; i < N; i++)
    {
       
        double x = std::rand() % (1);
        tep.set_vars({ {"x", &x} });
        /* Compile the expression with variables. */
        if (tep.compile(expression)) {
            double y = tep.evaluate();
        }
     
        cost_func.addObservation(x, y);
    }
 
    cost_func.solveByGaussNewton();
    return 0;
}


double simpson(double(*f)(double), double a, double b, int n)
{
    double s, dx, x;
    // if n is odd - add +1 interval to make it even
    if ((n / 2) * 2 != n) { n = n + 1; }
    s = 0.0;
    dx = (b - a) / static_cast<float>(n);
    for (int i = 2; i <= n - 1; i = i + 2)
    {
        x = a + static_cast<float>(i) * dx;
        s = s + 2.0 * f(x) + 4.0 * f(x + dx);
    }
    s = (s + f(a) + f(b) + 4.0 * f(a + dx)) * dx / 3.0;
    return s;
}





class LegendrePolynomial {
public:
    LegendrePolynomial(double lowerBound, double upperBound, size_t numberOfIterations)
        : mLowerBound(lowerBound), mUpperBound(upperBound), mNumberOfIterations(numberOfIterations), mWeight(numberOfIterations + 1), mRoot(numberOfIterations + 1) {
        calculateWeightAndRoot();
    }

    const std::vector<double>& getWeight() const {
        return mWeight;
    }

    const std::vector<double>& getRoot() const {
        return mRoot;
    }

private:
    const static double EPSILON;

    struct Result {
        double value;
        double derivative;

        Result() : value(0), derivative(0) {}
        Result(double val, double deriv) : value(val), derivative(deriv) {}
    };

    void calculateWeightAndRoot() {
        for (int step = 0; step <= mNumberOfIterations; step++) {
            double root = cos(M_PI * (step - 0.25) / (mNumberOfIterations + 0.5));
            Result result = calculatePolynomialValueAndDerivative(root);

            double newtonRaphsonRatio;
            do {
                newtonRaphsonRatio = result.value / result.derivative;
                root -= newtonRaphsonRatio;
                result = calculatePolynomialValueAndDerivative(root);
            } while (fabs(newtonRaphsonRatio) > EPSILON);

            mRoot[step] = root;
            mWeight[step] = 2.0 / ((1 - root * root) * result.derivative * result.derivative);
        }
    }

    Result calculatePolynomialValueAndDerivative(double x) {
        Result result(x, 0);

        double value_minus_1 = 1;
        const double f = 1 / (x * x - 1);
        for (int step = 2; step <= mNumberOfIterations; step++) {
            const double value = ((2 * step - 1) * x * result.value - (step - 1) * value_minus_1) / step;
            result.derivative = step * f * (x * value - result.value);

            value_minus_1 = result.value;
            result.value = value;
        }

        return result;
    }

    const double mLowerBound;
    const double mUpperBound;
    const int mNumberOfIterations;
    std::vector<double> mWeight;
    std::vector<double> mRoot;
};

const double LegendrePolynomial::EPSILON = 1e-15;

double gaussLegendreIntegral(double a, double b, int n, const std::function<double(double)>& f) {
    const LegendrePolynomial legendrePolynomial(a, b, n);
    const std::vector<double>& weight = legendrePolynomial.getWeight();
    const std::vector<double>& root = legendrePolynomial.getRoot();

    const double width = 0.5 * (b - a);
    const double mean = 0.5 * (a + b);

    double gaussLegendre = 0;
    for (int step = 1; step <= n; step++) {
        gaussLegendre += weight[step] * f(width * root[step] + mean);
    }

    return gaussLegendre * width;
}


double f(double x)
{
    return atan(sin(2 * x + 1));
}

// Takes lower bound, upper bound and inner point
double calculateExtremum(double a, double b, double x1, bool isMinimum)
{
    if ((b - a) < accuracy)	// End of searching
        return f((a + b) / 2);

    double x2;	// Establish points x1 and x2
    if (x1 > (a + b) / 2)	// x1 is actually x2
    {
        x2 = x1;
        x1 = b - (b - a) * goldenRatio;	// Calculate x1 based on the length of the interval
    }
    else
        x2 = a + (b - a) * goldenRatio; // Calculate x2 based on the length of the interval

    if (isMinimum)
    {
        if (f(x1) >= f(x2))
            return calculateExtremum(x1, b, x2, true);		// New interval is <x1, b> with x2 as inner point
        else
            return calculateExtremum(a, x2, x1, true);		// New interval is <a, x2> with x1 as inner point
    }
    else
    {
        if (f(x1) <= f(x2))
            return calculateExtremum(x1, b, x2, false);		// New interval is <x1, b> with x2 as inner point
        else
            return calculateExtremum(a, x2, x1, false);		// New interval is <a, x2> with x1 as inner point
    }

    return 0;
}

int golden()
{
    double potentialMinimum = calculateExtremum(lowerBound, upperBound, upperBound - ((upperBound - lowerBound) * goldenRatio), true);
    if (potentialMinimum < f(lowerBound) && potentialMinimum < f(upperBound))		// Check if is smaller than values on the boundaries
    {
        std::cout << "Minimum: " << potentialMinimum << std::endl;
        solve = potentialMinimum;
    }
    else
    {
        double maximum = calculateExtremum(lowerBound, upperBound, upperBound - ((upperBound - lowerBound) * goldenRatio), false);
        std::cout << "Maximum: " << maximum << std::endl;
        solve = maximum;
    }


    return 0;
}



std::string gen_form(int degree)
{
    std::string letters = "abcdefghijklmnopqrstuvwxyz";
    std::string poly = "";
    int i = 0;

    while (degree > 1)
    {
        poly = poly + letters[i] + "x^" + std::to_string(degree) + "+";
        i++;
        degree--;
    }
    poly = poly + letters[i] + "x+";
    i++;
    poly = poly + letters[i];

    return poly;
}


std::string polynomial(int degree, std::vector<int> coeffs)
{
    int i = 0;
    std::string poly = "";
    while (degree > 1)
    {
        if (coeffs[i] < 0 && i > 0)
        {
            poly.erase(poly.length() - 1, 1);
        }
        poly += std::to_string(coeffs[i]) + "x^" + std::to_string(degree) + "+";
        i++;
        degree--;
    }
    if (coeffs[i] < 0 && i > 0)
    {
        poly.erase(poly.length() - 1, 1);
    }
    poly += std::to_string(coeffs[i]) + "x" + "+";
    i++;
    if (coeffs[i] < 0 && i > 0)
    {
        poly.erase(poly.length() - 1, 1);
    }
    poly += std::to_string(coeffs[i]);
    return poly;
}


double p_of_x(double x, int degree, std::vector<int> coeffs)
{
    double value = 0;
    int i = 0;
    degree = degree + 1;
    while (degree--)
    {
        value = value + coeffs[i] * pow(x, degree);
        i++;
    }

    return value;
}


double derivative(double x, int degree, std::vector<int> coeffs)
{
    double value = 0;
    int i = 0;

    while (degree > 1)
    {
        value += degree * coeffs[i] * pow(x, degree - 1);
        i++;
        degree--;
    }
    value += coeffs[i];

    return value;
}


double* root(double x, int degree, std::vector<int> coeffs)
{
    double* y = new double[2];
    for (int i = 0; i < 200; i++)
    {
        double z = p_of_x(x, degree, coeffs);

        if (z > 0.000001 || z < -0.000001)
        {
            x = x - (p_of_x(x, degree, coeffs) / derivative(x, degree, coeffs));
        }
        else {
            y[0] = x;
            y[1] = 1;
            return y;
        }
    }
    y[0] = x;
    y[1] = 0;
    return y;
}

double Verify(double x, int degree, std::vector<int> coeffs)
{
    double value = p_of_x(x, degree, coeffs);
    value = (value * 1000000) / 100000;

    return value;
}


int newton_rapsody()
{
    int degree = 0;
    std::cout << "Enter degree of polynomial:";
    std::cin >> degree;
    std::string gen = gen_form(degree);
    std::cout << "General Form: " << gen << std::endl;

    int deg1 = degree + 1;
    std::vector<int> coeff;

    std::cout << "Enter coefficients(space separated): ";

    while (deg1--)
    {
        int temp;
        std::cin >> temp;
        coeff.push_back(temp);
    }
    std::cout << "Your polynomial is: " << polynomial(degree, coeff) << std::endl;
    double x;

    while (1)
    {
        std::cout << "Initial value: ";
        std::cin >> x;
        double* rt = root(x, degree, coeff);
        if (rt[1] == 1)
        {
            std::cout << "One of the roots(if multiple present): " << rt[0] << std::endl;
            solve  = rt[0];
        }
        else {
            std::cout << "No real roots found. The Loop stopped on value: " << rt[0] << std::endl;
            solve = rt[0];
        }
        std::cout << "Verififcation: P(found root or value loop stopped on) = " << Verify(rt[0], degree, coeff) << std::endl;
    }

}

double f2(double t, double x)
{
    const char* expression = name;
    double y;
    te_parser tep;
    tep.set_vars({ {"x", &x} });
    /* Compile the expression with variables. */
    if (tep.compile(expression)) {
        double y = tep.evaluate();
        int k = 0;
        return y;
    }


    double dx;
    dx = (-1.0) * x;
    //    dx = 0.5*sin(2*t) - x*cos(t);
    return dx;
}



double euler1d(double(*f)(double, double), double ti, double xi, double tf)
{
    double xf;
    xf = xi + f(ti, xi) * (tf - ti);
    return xf;
}


double euler1m(double(*f)(double, double), double ti, double xi, double tf)
{
    double xf;
    xf = xi + f(ti, xi) * (tf - ti);
    xf = xi + (f(ti, xi) + f(tf, xf)) * 0.5 * (tf - ti);
    return xf;
}


double rk4_1st(double(*f)(double, double), double ti, double xi, double tf)
{
    double xf;
    double h, k1, k2, k3, k4;

    h = tf - ti;

    k1 = h * f(ti, xi);
    k2 = h * f(ti + h / 2.0, xi + k1 / 2.0);
    k3 = h * f(ti + h / 2.0, xi + k2 / 2.0);
    k4 = h * f(ti + h, xi + k3);

    xf = xi + (k1 + 2.0 * (k2 + k3) + k4) / 6.0;
    return xf;
}



int euler_runge(int _key)
{
    double xi, ti, xf, tf, dt, tmax;
    int key = _key;
    const std::string method[3] = {"simple Euler","modified Euler","4th order Runge-Kutta"};

   

    /* initial information */
    key = 2;              // select a method (key = 0, 1, 2)   
    ti = 0.0;             // initial value for variable
    xi = 1.0;             // initial value for function
    dt = 0.1;             // step size for integration
    tmax = 2.0;          // inegrate from ti till tmax



    /* step 2: integration of ODE */
    while (ti <= tmax)
    {
        tf = ti + dt;
        if (key == 0) xf = euler1d(f2, ti, xi, tf);
        if (key == 1) xf = euler1m(f2, ti, xi, tf);
        if (key == 2) xf = rk4_1st(f2, ti, xi, tf);
        ti = tf;
        xi = xf;
        solve = xi;
    }
    return 0;
}























int main(int, char**)
{
    // Create application window
    //ImGui_ImplWin32_EnableDpiAwareness();
    WNDCLASSEXW wc = { sizeof(wc), CS_CLASSDC, WndProc, 0L, 0L, GetModuleHandle(NULL), NULL, NULL, NULL, NULL, L"ImGui Example", NULL };
    ::RegisterClassExW(&wc);
    HWND hwnd = ::CreateWindowW(wc.lpszClassName, L"Numerical Method", WS_OVERLAPPEDWINDOW, 100, 100, 1280, 800, NULL, NULL, wc.hInstance, NULL);

    // Initialize Direct3D
    if (!CreateDeviceD3D(hwnd))
    {
        CleanupDeviceD3D();
        ::UnregisterClassW(wc.lpszClassName, wc.hInstance);
        return 1;
    }

    // Show the window
    ::ShowWindow(hwnd, SW_SHOWDEFAULT);
    ::UpdateWindow(hwnd);

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsLight();

    // Setup Platform/Renderer backends
    ImGui_ImplWin32_Init(hwnd);
    ImGui_ImplDX9_Init(g_pd3dDevice);

    // Load Fonts
    // - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use ImGui::PushFont()/PopFont() to select them.
    // - AddFontFromFileTTF() will return the ImFont* so you can store it if you need to select the font among multiple.
    // - If the file cannot be loaded, the function will return NULL. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
    // - The fonts will be rasterized at a given size (w/ oversampling) and stored into a texture when calling ImFontAtlas::Build()/GetTexDataAsXXXX(), which ImGui_ImplXXXX_NewFrame below will call.
    // - Use '#define IMGUI_ENABLE_FREETYPE' in your imconfig file to use Freetype for higher quality font rendering.
    // - Read 'docs/FONTS.md' for more instructions and details.
    // - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
    //io.Fonts->AddFontDefault();
    //io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\segoeui.ttf", 18.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
    //ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f, NULL, io.Fonts->GetGlyphRangesJapanese());
    //IM_ASSERT(font != NULL);

    // Our state
    bool show_demo_window = false;
    bool show_another_window = false;
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    // Main loop
    bool done = false;
    while (!done)
    {
        // Poll and handle messages (inputs, window resize, etc.)
        // See the WndProc() function below for our to dispatch events to the Win32 backend.
        MSG msg;
        while (::PeekMessage(&msg, NULL, 0U, 0U, PM_REMOVE))
        {
            ::TranslateMessage(&msg);
            ::DispatchMessage(&msg);
            if (msg.message == WM_QUIT)
                done = true;
        }
        if (done)
            break;

        // Start the Dear ImGui frame
        ImGui_ImplDX9_NewFrame();
        ImGui_ImplWin32_NewFrame();
        ImGui::NewFrame();

        // 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
 

        // 2. Show a simple window that we create ourselves. We use a Begin/End pair to create a named window.
        {
            static float f = 0.0f;
            static int counter = 0;
     
            ImGui::Begin("Numerical Method Solver ");
      
           // Create a window called "Hello, world!" and append into it.
          //  

           // ImGui::LabelText();// Display some text (you can use a format strings too)
          //  ImGui::Checkbox("Demo Window", &show_demo_window);      // Edit bools storing our window open/close state
          //  ImGui::Checkbox("Another Window", &show_another_window);

       //     ImGui::SliderFloat("float", &f, 0.0f, 1.0f);            // Edit 1 float using a slider from 0.0f to 1.0f
       //     ImGui::ColorEdit3("clear color", (float*)&clear_color); // Edit 3 floats representing a color

            fn();
            ImGui::Text("Equation %s", name);

         //   ImGui::Text("Choose Solver.");
            
    //        double *coeff_arr = new double[argc - 4];	// Create a dynamically allocated coefficient array for the polynomial where coeff_arr[i] is equal to a_(n-i)
	//        for(int i = 1; i < argc - 3; i++)			// argv[0] is the executable filename so we store the command line arguments 1, 2, ..., argc - 4
	//	        coeff_arr[i-1] = atof(argv[i]);
	//
	//        double x_0 = atof(argv[argc - 3]);			// Remaining three command line arguments are x_0, x_1 and tolerance respectively.
	//        double x_1 = atof(argv[argc - 2]);			// We use the atof function to convert character arrays (strings) to floating point numbers.
	//        double tol = atof(argv[argc - 1]);
	//
	//        if(x_0 >= x_1) {							// x_0 < x_1 has to hold for the program to continue.
	//	        std::cerr << "[!ERROR!]: x_0 = " << x_0 << "  >=  x_1 = " << x_1 << std::endl;
	//	        return 1;				// Unsuccessful termination.
	//        }
	//
	//        if(tol <= 0.0) {							// Tolerance value has to be positive for the program to continue.
	//	        std::cerr << "[!ERROR!]: Tolerance value tol = " << tol << " has to be positive. " << std::endl;
	//	        return 1;				// Unsuccessful termination.
	//        }
	//
	//        // First number is the number of iterations, second number is the root for each pair.
	//
	//        // Print the calculated root and the number of iterations for the Bisection Method.
	//        std::pair<int, double> sol1 = bisectionMethod(coeff_arr, argc - 4, x_0, x_1, tol);
	//        std::cout << "For the Bisection Method: " << sol1.first << " number of iterations and " << sol1.second << " is the root." << std::endl;
            
          
            ImGui::Text("x = %.3f ", solve);
            ImGui::Text("write for simpson integration");
             char a[4] = "";
             char b[4] = "";
             char c[4] = "";
             if (ImGui::Button("Bisection"))
             {
                 bisection_caclulation();
             }
            if (ImGui::Button("Newton-Rapsody"))
            {
                newton_rapsody();
            }
            ImGui::InputText("Simpson a", a, IM_ARRAYSIZE(a));
            ImGui::InputText("Simpson b", b, IM_ARRAYSIZE(b));
            ImGui::InputText("Simpson n", c, IM_ARRAYSIZE(c));
            if (ImGui::Button("Simpson"))
            {
                double a_ = sscanf(a,"%d"  );
                double b_ = sscanf( b,"%d");
                double n_ = sscanf( c,"%n");
                simpson(fzero, a_, b_, n_);
            }
            ImGui::InputText("Legendrea a", a, IM_ARRAYSIZE(a));
            ImGui::InputText("Legendrea b", b, IM_ARRAYSIZE(b));
            ImGui::InputText("Legendrea n", c, IM_ARRAYSIZE(c));
            if (ImGui::Button("GaussLegendrea"))
            {
                double a_ = sscanf(a, "%d");
                double b_ = sscanf(b, "%d");
                double n_ = sscanf(c, "%n");
                gaussLegendreIntegral( a_, b_, n_, fzero);
            }
            if (ImGui::Button("GoldenDivide"))
            {
                golden();
            }
            if (ImGui::Button("NewtonOptimalization"))
            {
                newton();
            }
            if (ImGui::Button("Euler"))
            {
                euler_runge(1);
            }
            if (ImGui::Button("RungegoKutty"))
            {
                euler_runge(2);
            }
          

            ImGui::End();
        }

        // 3. Show another simple window.
        if (show_another_window)
        {
            ImGui::Begin("Another Window", &show_another_window);   // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
            ImGui::Text("Hello from another window!");
            if (ImGui::Button("Close Me"))
                show_another_window = false;
            ImGui::End();
        }

        // Rendering
        ImGui::EndFrame();
        g_pd3dDevice->SetRenderState(D3DRS_ZENABLE, FALSE);
        g_pd3dDevice->SetRenderState(D3DRS_ALPHABLENDENABLE, FALSE);
        g_pd3dDevice->SetRenderState(D3DRS_SCISSORTESTENABLE, FALSE);
        D3DCOLOR clear_col_dx = D3DCOLOR_RGBA((int)(clear_color.x*clear_color.w*255.0f), (int)(clear_color.y*clear_color.w*255.0f), (int)(clear_color.z*clear_color.w*255.0f), (int)(clear_color.w*255.0f));
        g_pd3dDevice->Clear(0, NULL, D3DCLEAR_TARGET | D3DCLEAR_ZBUFFER, clear_col_dx, 1.0f, 0);
        if (g_pd3dDevice->BeginScene() >= 0)
        {
            ImGui::Render();
            ImGui_ImplDX9_RenderDrawData(ImGui::GetDrawData());
            g_pd3dDevice->EndScene();
        }
        HRESULT result = g_pd3dDevice->Present(NULL, NULL, NULL, NULL);

        // Handle loss of D3D9 device
        if (result == D3DERR_DEVICELOST && g_pd3dDevice->TestCooperativeLevel() == D3DERR_DEVICENOTRESET)
            ResetDevice();
    }

    ImGui_ImplDX9_Shutdown();
    ImGui_ImplWin32_Shutdown();
    ImGui::DestroyContext();

    CleanupDeviceD3D();
    ::DestroyWindow(hwnd);
    ::UnregisterClassW(wc.lpszClassName, wc.hInstance);

    return 0;
}

// Helper functions

bool CreateDeviceD3D(HWND hWnd)
{
    if ((g_pD3D = Direct3DCreate9(D3D_SDK_VERSION)) == NULL)
        return false;

    // Create the D3DDevice
    ZeroMemory(&g_d3dpp, sizeof(g_d3dpp));
    g_d3dpp.Windowed = TRUE;
    g_d3dpp.SwapEffect = D3DSWAPEFFECT_DISCARD;
    g_d3dpp.BackBufferFormat = D3DFMT_UNKNOWN; // Need to use an explicit format with alpha if needing per-pixel alpha composition.
    g_d3dpp.EnableAutoDepthStencil = TRUE;
    g_d3dpp.AutoDepthStencilFormat = D3DFMT_D16;
    g_d3dpp.PresentationInterval = D3DPRESENT_INTERVAL_ONE;           // Present with vsync
    //g_d3dpp.PresentationInterval = D3DPRESENT_INTERVAL_IMMEDIATE;   // Present without vsync, maximum unthrottled framerate
    if (g_pD3D->CreateDevice(D3DADAPTER_DEFAULT, D3DDEVTYPE_HAL, hWnd, D3DCREATE_HARDWARE_VERTEXPROCESSING, &g_d3dpp, &g_pd3dDevice) < 0)
        return false;

    return true;
}

void CleanupDeviceD3D()
{
    if (g_pd3dDevice) { g_pd3dDevice->Release(); g_pd3dDevice = NULL; }
    if (g_pD3D) { g_pD3D->Release(); g_pD3D = NULL; }
}

void ResetDevice()
{
    ImGui_ImplDX9_InvalidateDeviceObjects();
    HRESULT hr = g_pd3dDevice->Reset(&g_d3dpp);
    if (hr == D3DERR_INVALIDCALL)
        IM_ASSERT(0);
    ImGui_ImplDX9_CreateDeviceObjects();
}

// Forward declare message handler from imgui_impl_win32.cpp
extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

// Win32 message handler
// You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
// - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application, or clear/overwrite your copy of the mouse data.
// - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application, or clear/overwrite your copy of the keyboard data.
// Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
LRESULT WINAPI WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
    if (ImGui_ImplWin32_WndProcHandler(hWnd, msg, wParam, lParam))
        return true;

    switch (msg)
    {
    case WM_SIZE:
        if (g_pd3dDevice != NULL && wParam != SIZE_MINIMIZED)
        {
            g_d3dpp.BackBufferWidth = LOWORD(lParam);
            g_d3dpp.BackBufferHeight = HIWORD(lParam);
            ResetDevice();
        }
        return 0;
    case WM_SYSCOMMAND:
        if ((wParam & 0xfff0) == SC_KEYMENU) // Disable ALT application menu
            return 0;
        break;
    case WM_DESTROY:
        ::PostQuitMessage(0);
        return 0;
    }
    return ::DefWindowProc(hWnd, msg, wParam, lParam);
}
