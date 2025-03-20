// MACHINE PROBLEM (1)
function acp_compute() {
    try {
        
        if (document.getElementById("acp_variable").value == "" || document.getElementById("acp_point").value == "" || document.getElementById("acp_degree").value == "" || document.getElementById("acp_decimal").value== "") {
            window.alert("Input incomplete. Please fill up the required fields.")
        }
        else {
            variable    = parseFloat(document.getElementById("acp_variable").value);
            point       = parseFloat(document.getElementById("acp_point").value);
            degree      = parseFloat(document.getElementById("acp_degree").value);
            decimal     = parseFloat(document.getElementById("acp_decimal").value);

            //Formula
            dx_cos   = -Math.sin(point);
            dx_nsin  = -Math.cos(point);
            dx_ncos  =  Math.sin(point);
            dx_sin   =  Math.cos(point);

            res_fm = [];
            deg_fm = degree;

            while(deg_fm > 0) {
                if(deg_fm > 0 ){
                    res_fm.push(dx_cos);
                    deg_fm -= 1;
                }
                if(deg_fm > 0 ){
                    res_fm.push(dx_nsin);
                    deg_fm -= 1;
                }
                if(deg_fm > 0 ){
                    res_fm.push(dx_ncos);
                    deg_fm -= 1;
                }
                if(deg_fm > 0 ){
                    res_fm.push(dx_sin);
                    deg_fm -= 1;
                }
            }


            // Final Calculation
            final_ans = Math.cos(point);
            num_rise  = 0;

            for( i = 0; i < degree; i++ ) {
                num_rise += 1;
                if (point == 0) {
                    final_ans += ((res_fm[i]) * variable ** num_rise) / factorial(i+1);
                }
                else {
                    final_ans += ( (1 / factorial(i+1)) * (variable - point) ** num_rise) * (res_fm[i])
                }
            }

            // Results
            document.getElementById("result_true_value").innerHTML          = Math.cos(variable);
            document.getElementById("result_approximate_vlaue").innerHTML   = final_ans;
            document.getElementById("result_chopping").innerHTML            = truncate(final_ans, decimal);
            document.getElementById("result_rounding").innerHTML            = final_ans.toFixed(decimal);

            // Errors
            document.getElementById("ate_av").innerHTML = Math.abs(Math.cos(variable) - final_ans);
            document.getElementById("ate_c").innerHTML  = Math.abs(Math.cos(variable) - truncate(final_ans, decimal));
            document.getElementById("ate_r").innerHTML  = Math.abs(Math.cos(variable) - final_ans.toFixed(decimal));

            document.getElementById("pe_av").innerHTML = (Math.abs(Math.cos(variable) - final_ans) / Math.cos(variable)) * 100 + "%";
            document.getElementById("pe_c").innerHTML  = (Math.abs(Math.cos(variable) - truncate(final_ans, decimal)) / Math.cos(variable)) * 100 + "%";
            document.getElementById("pe_r").innerHTML  = (Math.abs(Math.cos(variable) - final_ans.toFixed(decimal)) / Math.cos(variable)) * 100 + "%";
        }
    }
    catch(e) {
        window.alert("There was an error in your input field. Please try again.")
    }
}
function factorial(num) {
    if (num < 0) 
          return -1;
    else if (num == 0) 
        return 1;
    else {
        return (num * factorial(num - 1));
    }
}
function truncate(number, digits) {
    stepper = 10.0 ** digits
    return Math.trunc(stepper * number) / stepper
}
function acp_clear() {
    document.getElementById("acp_variable").value   = "";
    document.getElementById("acp_point").value      = "";
    document.getElementById("acp_degree").value     = "";
    document.getElementById("acp_decimal").value    = "";

    // Results
    document.getElementById("result_true_value").innerHTML          = "---";
    document.getElementById("result_approximate_vlaue").innerHTML   = "---";
    document.getElementById("result_chopping").innerHTML            = "---";
    document.getElementById("result_rounding").innerHTML            = "---";

    // Errors
    document.getElementById("ate_av").innerHTML = "---";
    document.getElementById("ate_c").innerHTML  = "---";
    document.getElementById("ate_r").innerHTML  = "---";

    document.getElementById("pe_av").innerHTML = "---";
    document.getElementById("pe_c").innerHTML  = "---";
    document.getElementById("pe_r").innerHTML  = "---";
}

// MACHINE PROBLEM (2)
function solve_bisection() {
    //Try
    try {
      //if formula is null
      if(document.getElementById("formula").value == "") {
        window.confirm("Incomplete: Enter formula (f(x))");
      }
      //if interval a is null
      else if(document.getElementById("value_a").value == "") {
        window.confirm("Incomplete: Enter interval [a]");
      }
      //if interval b is null
      else if(document.getElementById("value_b").value == "") {
        window.confirm("Incomplete: Enter interval [b]");
      }
      //if error is null
      else if(document.getElementById("error").value == "") {
        window.confirm("Incomplete: Enter error (ε>0)");
      }
      //if complete = evaluate
      else {
        //Get value of interval [a] and [b]; with formula for c
        a = parseFloat(document.getElementById("value_a").value);
        b = parseFloat(document.getElementById("value_b").value);
        c = (a+b) / 2;
        //Get value of Error | Epsilon
        epsilon = parseFloat(document.getElementById("error").value);
        //Get formula 
        formula_raw = document.getElementById("formula").value;
        //CONVERSION (symbols of polynomials, trigonometric, and exponential)
        formula_sign = formula_raw.replace(/pi/g,"Math.PI").replace(/0x/g,"0").replace(/1x/g,"x").replace(/2x/g,"2*x").replace(/3x/g,"3*x").replace(/4x/g,"4*x").replace(/5x/g,"5*x").replace(/6x/g,"6*x").replace(/7x/g,"7*x").replace(/8x/g,"8*x").replace(/9x/g,"9*x");
        formula_trig = formula_sign.replace(/sin/g,"Math.sin").replace(/cos/g,"Math.cos").replace(/tan/g,"Math.tan").replace(/csc/g,"1/Math.sin").replace(/sec/g,"1/Math.cos").replace(/cot/g,"1/Math.tan");
        formula_exp  = formula_trig.replaceAll("^","**").replace(/e/g,"Math.E");
        formula_html = formula_exp;
        //Convert letters
        formula_a = formula_html.replace(/x/g,"a");
        formula_b = formula_html.replace(/x/g,"b");
        formula_c = formula_html.replace(/x/g,"c");
        //Convert into expression
        expression_a = eval(formula_a);
        expression_b = eval(formula_b);
        expression_c = eval(formula_c);
        //Array containers for values
        n   = [];
        a_a = [];
        b_b = [];
        fa  = [];
        fb  = [];
        c_c = [];
        fc  = [];
        final = [];
        //Push first iteration
        //n
        n.push(1)
        //an
        a_a.push(a);
        //bn
        b_b.push(b);
        //cn
        c_c.push(c);
        //f(a)
        fa.push(expression_a);
        //f(b)
        fb.push(expression_b);
        //f(c)
        fc.push(expression_c);
        //|bn - an| / 2
        final.push(Math.abs(b - a)/2)
        // ASSUMPTION < 0
        if((fa[0]*fb[0])<0){
          //START ITERATION -------------------------------------
          //Assumptions if Bisection is possible
          document.getElementById("a1").innerHTML = "f(x) is continous on [" + a + "," + b + "]\n and f(" + a + ") * f(" + b + ") = " + "(" + fa[0] + ")(" + fb[0] + ") = " + (fa[0]*fb[0]) + " < 0";
          document.getElementById("a2").innerHTML = "> Assumptions are satisfied";
          document.getElementById("a3").innerHTML = "> Bisection can be used";
          whl = 1;
          i = 0
          while(whl == 1){
            //Push n values
            //n
            n.push(i+2);
            //if fa(n) * fc(n) < 0
            if((fa[i]*fc[i])<0){
              b = c;
              //push c(n-1) -> b(n)
              b_b.push(c_c[i]);
            }
            else{
              //push b(n-1) -> b(n)
              b_b.push(b_b[i]);
            }
            //if fa(n) * fc(n) > 0
            if((fa[i]*fc[i])>0){
              a = c
              //push c(n-1) -> a(n)
              a_a.push(c_c[i]);
            }
            else{
              //push a(n-1) -> a(n)
              a_a.push(a_a[i]);
            }
            //push c(n)
            c_c.push((a_a[i+1] + b_b[i+1]) / 2); 
            //push f(c)
            c = (a+b)/2;
            // push evaluate f(a)
            fa.push(eval(formula_a));
            // push evaluate f(b)
            fb.push(eval(formula_b));
            // push evaluate f(c)
            fc.push(eval(formula_c))
            // push |bn - an| / 2
            final.push(Math.abs(b - a)/2)
            // if (|bn - an| / 2) < error stop iteration
            if(epsilon >= final[i+1]){
              whl = 0;
            }
            // Max iteration if loop is infinite
            else if(n.length >= 1000) {
              whl = 0;
            }
            i++;
          }
          //------------------------------------------------
          //Table of values
          const tbl = [n, a_a, b_b, fa, fb, c_c, fc, final];
          //Execute print table
          var table = document.getElementsByTagName('table')[4];
          //loop table values from array
          for(var i = 0; i <n.length;i++) {
            var newRow = table.insertRow(table.rows.length); 
            var cel1 = newRow.insertCell(0);
            var cel2 = newRow.insertCell(1);
            var cel3 = newRow.insertCell(2);
            var cel4 = newRow.insertCell(3);
            var cel5 = newRow.insertCell(4);
            var cel6 = newRow.insertCell(5);
            var cel7 = newRow.insertCell(6);
            var cel8 = newRow.insertCell(7);
            // input the values from the array to the table
            cel1.innerHTML = n[i];
            cel2.innerHTML = a_a[i];
            cel3.innerHTML = b_b[i];
            cel4.innerHTML = fa[i];
            cel5.innerHTML = fb[i];
            cel6.innerHTML = c_c[i];
            cel7.innerHTML = fc[i];
            cel8.innerHTML = final[i];
          }
        }
        // ASSUMPTION > 0
        else{
          // Assumptions if Bisection is not possible
          document.getElementById("a1").innerHTML = "f(x) is continous on [" + a + "," + b + "]\n and f(" + a + ") * f(" + b + ") = " + "(" + fa[0] + ")(" + fb[0] + ") = " + (fa[0]*fb[0]) + " > 0";
          document.getElementById("a2").innerHTML = "> Assumptions are NOT satisfied";
          document.getElementById("a3").innerHTML = "> Bisection CANNOT be used";
        }
      }
    } 
    //Catch error exemption
    catch (e) {
        window.confirm("There is an error in your input. Please try again.")
    }
}
function solve_newton() {
    //Try
    try {
      //if formula is null
      if(document.getElementById("formula_n").value == "") {
        window.confirm("Incomplete: Enter formula (f(x))");
      }
      //if initial point x0 is null
      else if(document.getElementById("initial_point").value == "") {
        window.confirm("Incomplete: Enter initial point (x0)");
      }
      //if error is null
      else if(document.getElementById("error_n").value == "") {
        window.confirm("Incomplete: Enter error (ε>0)");
      }
      //if complete = evaluate
      else {
        //Get Initial point Xn
        initial_point = parseFloat(document.getElementById("initial_point").value);
        x = initial_point;
        //Get value of Error | Epsilon
        epsilon = parseFloat(document.getElementById("error_n").value);
          //Get formula 
          formula_raw = document.getElementById("formula_n").value;
          formula_raw_derivative = String(math.derivative(formula_raw, 'x'));
          //CONVERSION (symbols of polynomials, trigonometric, and exponential)
          formula_sign = formula_raw.replace(/pi/g,"Math.PI").replace(/0x/g,"0").replace(/1x/g,"x").replace(/2x/g,"2*x").replace(/3x/g,"3*x").replace(/4x/g,"4*x").replace(/5x/g,"5*x").replace(/6x/g,"6*x").replace(/7x/g,"7*x").replace(/8x/g,"8*x").replace(/9x/g,"9*x");
          formula_trig = formula_sign.replace(/sin/g,"Math.sin").replace(/cos/g,"Math.cos").replace(/tan/g,"Math.tan").replace(/csc/g,"1/Math.sin").replace(/sec/g,"1/Math.cos").replace(/cot/g,"1/Math.tan");
          formula_exp  = formula_trig.replaceAll("^","**").replace(/e/g,"Math.E");
          formula_html = formula_exp;
          //CONVERSION formula - derivative
          formula_sign_d = formula_raw_derivative.replace(/pi/g,"Math.PI").replace(/0x/g,"0").replace(/1x/g,"x").replace(/2x/g,"2*x").replace(/3x/g,"3*x").replace(/4x/g,"4*x").replace(/5x/g,"5*x").replace(/6x/g,"6*x").replace(/7x/g,"7*x").replace(/8x/g,"8*x").replace(/9x/g,"9*x");
          formula_trig_d = formula_sign_d.replace(/sin/g,"Math.sin").replace(/cos/g,"Math.cos").replace(/tan/g,"Math.tan").replace(/csc/g,"1/Math.sin").replace(/sec/g,"1/Math.cos").replace(/cot/g,"1/Math.tan");
          formula_exp_d  = formula_trig_d.replaceAll("^","**").replace(/e/g,"Math.E");
          formula_html_d = formula_exp_d;
          //Convert into expression
          expression_xn = eval(formula_html);
          expression_fxn = eval(formula_html_d);
        //array containers
        n   = [];
        xn  = [];
        fxn = [];
        fxn_d = [];
        xn_1  = [];
        final = [];
        //First Iteration
        //push n 
        n.push(1);
        //push initial point x0
        xn.push(initial_point);
        //push f(xn)
        fxn.push(expression_xn);
        //push f'(xn)
        fxn_d.push(expression_fxn);
        //psuh xn+1
        xn_1.push(xn - (fxn[0] / fxn_d[0]));
        //push |xn+1 - xn|
        final.push(Math.abs(xn_1[0] - xn[0]));
        //Table of values
        table = [n, xn, fxn, fxn_d, xn_1, final];
        //if f(xn) != 0
        if(fxn_d[0] != 0){
          //Assumptions if Newton-Raphson is possible
          document.getElementById("a1_n").innerHTML = "f(x) is continous and the first derivative is known. f'(x) = " + math.derivative(formula_raw, 'x');
          document.getElementById("a2_n").innerHTML = "An inital guess x0 such that f'(x) ≠ 0 is given. f'(x) = " + fxn_d[0];
          document.getElementById("a3_n").innerHTML = "> Assumptins are satisfied";
          document.getElementById("a4_n").innerHTML = "> Newton can be used";
          i = 0;
          whl = 1
          //loop iteration 
          while(whl==1) {
            // push n
            n.push(i+2);
            // Change
            x = xn_1[i];
            // push xn
            xn.push(x);
            // push f(xn)
            fxn.push(eval(formula_html));
            // push f'(xn)
            fxn_d.push(eval(formula_html_d));
            // push (xn+1)
            xn_1.push(xn[i+1] - (fxn[i+1] / fxn_d[i+1]));
            // push (|xn+1 - xn|)
            final.push(Math.abs(xn_1[i+1] - xn[i+1]));
            //if (|xn+1 - xn|) <= error stop iteration
            if(epsilon >= final[i+1]){
              whl = 0;
            }
            // Max iteration if loop is infinite
            else if(n.length >= 1000) {
              whl = 0;
            }
            i++;
          }
          // execute table
          var table = document.getElementsByTagName('table')[5];
          // loop the values of the table
          for(var i = 0; i <n.length;i++) {
            var newRow = table.insertRow(table.rows.length); 
            // add rows to the table
            var cel1 = newRow.insertCell(0);
            var cel2 = newRow.insertCell(1);
            var cel3 = newRow.insertCell(2);
            var cel4 = newRow.insertCell(3);
            var cel5 = newRow.insertCell(4);
            var cel6 = newRow.insertCell(5);
            // input the values from the array to the table
            cel1.innerHTML = n[i];
            cel2.innerHTML = xn[i];
            cel3.innerHTML = fxn[i];
            cel4.innerHTML = fxn_d[i];
            cel5.innerHTML = xn_1[i];
            cel6.innerHTML = final[i];
          }
        }
        else{
          // Assumptions if Newton-Raphson is not possible
          document.getElementById("a1_n").innerHTML = "An inital guess x0 such that f'(x) ≠ 0 is given. f'(x) = " + fxn_d[0];
          document.getElementById("a2_n").innerHTML = "> Assumptins are NOT satisfied";
          document.getElementById("a3_n").innerHTML = "> Newton CANNOT be used";
        }
      }
    }
    //Catch error excemption
    catch (e) {
      window.confirm("There is an error in your input. Please try again.")
    }
}

// MACHINE PROBLEM (3)
// store points entered
arr_x = [];
arr_y = [];
// Dollar to PHP exchange rate
// LINK: https://tradingeconomics.com/philippines/currency
// Days
arr_x_pd = [1,6,11,12,14,15,19,20,21,22,25,27,28,29,32,34,39,41,42,43];
// Value of dollars
arr_y_pd = [51.52,51.27,52.1,52.02,52.18,52.15,52.46,52.35,52.43,52.33,52.38,52.15,52.37,52.31,52.54,52.39,52.589,52.17,52.38,52.37];
// Add function
function add() {
    try {
        // Retrieve the values form the input field
        x = document.getElementById("value_x").value
        y = document.getElementById("value_y").value
        // if the input field x and y are not values
        if(isNaN(x) | isNaN(y)) {
            window.confirm("There is an error in your input. Please try again.");
            // clear input 
            document.getElementById("value_x").value = "";
            document.getElementById("value_y").value = "";
        }
        // else if it is empty
        else if(x == "" | y == "") {
            window.confirm("Please input points. Try again.");
        }
        // if requirements are satisfied then initialize
        else {
            // Push the values to the array
            arr_x.push(x);
            arr_y.push(y);
            // Display all the values
            display = []
            for(i = 0; i < arr_x.length; i++) {
                display.push("(" + arr_x[i] + "," + arr_y[i] + ")")
            }
            document.getElementById("points").innerHTML = display;
            // clear input 
            document.getElementById("value_x").value = "";
            document.getElementById("value_y").value = "";
        }
    }
    catch (e) {
        window.confirm("There is an error in your input. Please try again.");
    }
}
// Solve functions
function solve() {
    try {
        // Retrieve value of x
        x = document.getElementById("value_for_x").value;
        // if the input field x is not a values
        if (isNaN(x)) {
            window.confirm("There is an error in your input. Please try again.");
            // clear input 
            document.getElementById("value_for_x").value = "";
        }
        // else if it is empty
        else if (x == ""){
            window.confirm("Please input a value for x. Try again.");
        }
        // if requirements are satisfied then initialize
        else {
            // Lagrange interpolation
            store_values = [];
            for (i = 0; i < arr_x.length; i++) {
                for (j = 0; j < arr_y.length; j++) {
                    if (j != i) {
                        store_values.push( ( (x - arr_x[j]) / (arr_x[i] - arr_x[j]) ) );
                    }
                }
            }
            // Group values 
            product_values = [];
            product = 1;
            index_store_values = 0;
            for (i = 0; i < arr_x.length; i++) {
                for (j = 0; j < arr_y.length-1; j++) {
                    product *= store_values[index_store_values];
                    index_store_values ++;
                }
                product_values.push(product);
                product = 1;
            }
            // Sum all the values
            sum = 0;
            for (i = 0; i < product_values.length; i++) {
                sum += product_values[i] * arr_y[i];
            }
            // For Value testing visualization
            // document.getElementById("answer1").innerHTML = store_values;
            // document.getElementById("answer2").innerHTML = product_values;
            document.getElementById("answer").innerHTML = "f(" + x + ") = " + sum;
            // clear input 
            document.getElementById("value_for_x").value = "";
            document.getElementById("sub_answer").value = product_values;
        }
    }
    catch (e) {
        
    }
}
function solve1() {
    try {
        // Retrieve value of x
        x = document.getElementById("value_for_x_pd").value;
        // if the input field x is not a values
        if (isNaN(x)) {
            window.confirm("There is an error in your input. Please try again.");
            // clear input 
            document.getElementById("value_for_x").value = "";
        }
        // else if it is empty
        else if (x == ""){
            window.confirm("Please input a value for x. Try again.");
        }
        // if requirements are satisfied then initialize
        else {
            // Lagrange interpolation
            store_values = [];
            for (i = 0; i < arr_x_pd.length; i++) {
                for (j = 0; j < arr_y_pd.length; j++) {
                    if (j != i) {
                        store_values.push( ( (x - arr_x_pd[j]) / (arr_x_pd[i] - arr_x_pd[j]) ) );
                    }
                }
            }
            // Group values 
            product_values = [];
            product = 1;
            index_store_values = 0;
            for (i = 0; i < arr_x_pd.length; i++) {
                for (j = 0; j < arr_y_pd.length-1; j++) {
                    product *= store_values[index_store_values];
                    index_store_values ++;
                }
                product_values.push(product);
                product = 1;
            }
            // Sum all the values
            sum = 0;
            for (i = 0; i < product_values.length; i++) {
                sum += product_values[i] * arr_y_pd[i];
            } 
            // For Value testing visualization
            // document.getElementById("answer1").innerHTML = store_values;
            // document.getElementById("answer2").innerHTML = arr_x_pd.length;
            document.getElementById("answer_pd").innerHTML = "f(" + x + ") = " + sum;
            // clear input 
            document.getElementById("sub_answer").value = product_values;
        }
    }
    catch (e) {
        
    }
}
// Clear function
function clear_all() {
    // Clear all values of array
    arr_x = [];
    arr_y = [];
    // clear display
    document.getElementById("points").innerHTML = "";
    document.getElementById("answer").innerHTML = ""; 
    // clear input 
    document.getElementById("value_x").value = "";
    document.getElementById("value_y").value = "";
    document.getElementById("value_for_x").value = "";
    // For Value testing visualization
    // document.getElementById("answer1").innerHTML = ""; 
    // document.getElementById("answer2").innerHTML = ""; 
}
function clear_all1() {
    // clear display
    document.getElementById("answer_pd").innerHTML = ""; 
    // clear input 
    document.getElementById("value_for_x_pd").value = "";
    // For Value testing visualization
    // document.getElementById("answer1").innerHTML = ""; 
    // document.getElementById("answer2").innerHTML = ""; 
}

// MACHINE PROBLEM (4)
function compute_trapezoidal_simpson() {
    func        = document.getElementById("value_function").value
    a           = parseFloat(document.getElementById("value_a_mp4").value);
    b           = parseFloat(document.getElementById("value_b_mp4").value)
    n           = parseFloat(document.getElementById("subinterval").value)
    delta_x     = parseFloat((b - a) / n);

    try {
        if ( document.getElementById("value_function").value == "" || document.getElementById("value_a_mp4").value == "" || document.getElementById("value_b_mp4").value == "" || document.getElementById("subinterval").value == "") {
            window.alert("Input incomplete. Please fill up all given inputs")
        }
        else {
            if (b > a) {
                //CONVERSION (symbols of polynomials, trigonometric, and exponential)
                formula_sign = func.replace(/pi/g,"Math.PI").replace(/sqrt/g,"Math.sqrt").replace(/0x/g,"0").replace(/1x/g,"x").replace(/2x/g,"2*x").replace(/3x/g,"3*x").replace(/4x/g,"4*x").replace(/5x/g,"5*x").replace(/6x/g,"6*x").replace(/7x/g,"7*x").replace(/8x/g,"8*x").replace(/9x/g,"9*x");
                formula_trig = formula_sign.replace(/sin/g,"Math.sin").replace(/cos/g,"Math.cos").replace(/tan/g,"Math.tan").replace(/csc/g,"1/Math.sin").replace(/sec/g,"1/Math.cos").replace(/cot/g,"1/Math.tan");
                formula_exp  = formula_trig.replaceAll("^","**").replace(/e/g,"Math.E");
                formula_html = formula_exp;
    
                // initial value xi
                store_values_xn = [];
                for (i = 0; i < n; i++) {
                    store_values_xn.push(a + (delta_x * (i+1))); // xi
                }
    
                // function // store the values
                store_values_xn_function = [];
                x = a;
                store_values_xn_function.push(eval(formula_html)); // for x0
                for (i = 0; i < n; i++) {
                    x = store_values_xn[i];
                    store_values_xn_function.push(eval(formula_html)); // for xn 
                }
    
                // multiply by 2 - trapezoidal
                store_values_xn_sum_trapezoidal = 0;
                for (i = 0; i < store_values_xn_function.length; i++) {
                    if(i == 0 || i == store_values_xn_function.length - 1) {
                        store_values_xn_sum_trapezoidal += store_values_xn_function[i];
                    }
                    else {
                        store_values_xn_sum_trapezoidal += store_values_xn_function[i] * 2;
                    }
                }
    
                // multiply by 2 and 4 - simpson
                store_values_xn_sum_simpson = 0
                for (i = 0; i < store_values_xn_function.length; i++) {
                    if(i == 0 || i == store_values_xn_function.length - 1) {
                        store_values_xn_sum_simpson += store_values_xn_function[i];
                    }
                    else {
                        if(i % 2 == 0) {
                            store_values_xn_sum_simpson += store_values_xn_function[i] * 2;
                        }
                        else {
                            store_values_xn_sum_simpson += store_values_xn_function[i] * 4;
                        }
                    }
                }
    
                // Trapezoidal
                trapezoidal = (delta_x / 2) * store_values_xn_sum_trapezoidal;
                // Simpson
                simpson = (delta_x / 3) * store_values_xn_sum_simpson;
    
                // Print 
                document.getElementById("result_trapezoidal").innerHTML = trapezoidal;
                if (n % 2 == 0) {
                    document.getElementById("result_simpson").innerHTML = simpson;
                }
                else {
                    document.getElementById("result_simpson").innerHTML = "subintevral (n) should be an even number"
                }
            }
            else {
                window.alert("The lower limit (a) should be less than the upper limit (b).")
            }
        }
    }
    catch (e) {
        window.alert("There is an error in your input. Please try again.")
    }

    // x^2-2x+2
    // 2x^2-4x+1
}
function clear_trapezoidal_simpson() {
    document.getElementById("value_function").value         = "";
    document.getElementById("value_a_mp4").value                = "";
    document.getElementById("value_b_mp4").value                = "";
    document.getElementById("subinterval").value            = "";
    document.getElementById("result_trapezoidal").innerHTML = "";
    document.getElementById("result_simpson").innerHTML     = "";
}

// functions
var MP1   = document.getElementById("MP1");
var MP2_1 = document.getElementById("MP2_1");
var MP2_2 = document.getElementById("MP2_2");
var MP3_1 = document.getElementById("MP3_1");
var MP3_2 = document.getElementById("MP3_2");
var MP4   = document.getElementById("MP4")

//Reset design
function MP_btn() {
    document.getElementById("MP1_btn").style.backgroundColor = "";
    document.getElementById("MP1_btn").style.color = "";

    document.getElementById("MP2_1_btn").style.backgroundColor = "";
    document.getElementById("MP2_1_btn").style.color = "";

    document.getElementById("MP2_2_btn").style.backgroundColor = "";
    document.getElementById("MP2_2_btn").style.color = "";

    document.getElementById("MP3_1_btn").style.backgroundColor = "";
    document.getElementById("MP3_1_btn").style.color = "";

    document.getElementById("MP3_2_btn").style.backgroundColor = "";
    document.getElementById("MP3_2_btn").style.color = "";

    document.getElementById("MP4_btn").style.backgroundColor = "";
    document.getElementById("MP4_btn").style.color = "";
}
//Show the MP's
function MP1_btn() {
    MP_btn();
    if (MP1.style.display === "none") {
        MP1.style.display = "block";
        MP2_1.style.display = "none";
        MP2_2.style.display = "none";
        MP3_1.style.display = "none";
        MP3_2.style.display = "none";
        MP4.style.display = "none";
        document.getElementById("MP1_btn").style.backgroundColor = "#4CAF50";
        document.getElementById("MP1_btn").style.color = "white";
    } else {
        MP1.style.display = "none";
        document.getElementById("MP1_btn").style.backgroundColor = "";
        document.getElementById("MP1_btn").style.color = "";
    }
}
function MP2_1_btn() {
    MP_btn();
    if (MP2_1.style.display === "none") {
        MP1.style.display = "none";
        MP2_1.style.display = "block";
        MP2_2.style.display = "none";
        MP3_1.style.display = "none";
        MP3_2.style.display = "none";
        MP4.style.display = "none";
        document.getElementById("MP2_1_btn").style.backgroundColor = "#4CAF50";
        document.getElementById("MP2_1_btn").style.color = "white";
    } else {
        MP2_1.style.display = "none";
        document.getElementById("MP2_1_btn").style.backgroundColor = "";
        document.getElementById("MP2_1_btn").style.color = "";
    }
}
function MP2_2_btn() {
    MP_btn();
    if (MP2_2.style.display === "none") {
        MP1.style.display = "none";
        MP2_1.style.display = "none";
        MP2_2.style.display = "block";
        MP3_1.style.display = "none";
        MP3_2.style.display = "none";
        MP4.style.display = "none";
        document.getElementById("MP2_2_btn").style.backgroundColor = "#4CAF50";
        document.getElementById("MP2_2_btn").style.color = "white";
    } else {
        MP2_2.style.display = "none";
        document.getElementById("MP2_2_btn").style.backgroundColor = "";
        document.getElementById("MP2_2_btn").style.color = "";
    }
}
function MP3_1_btn() {
    MP_btn();
    if (MP3_1.style.display === "none") {
        MP1.style.display = "none";
        MP2_1.style.display = "none";
        MP2_2.style.display = "none";
        MP3_1.style.display = "block";
        MP3_2.style.display = "none";
        MP4.style.display = "none";
        document.getElementById("MP3_1_btn").style.backgroundColor = "#4CAF50";
        document.getElementById("MP3_1_btn").style.color = "white";
    } else {
        MP3_1.style.display = "none";
        document.getElementById("MP3_1_btn").style.backgroundColor = "";
        document.getElementById("MP3_1_btn").style.color = "";
    }
}
function MP3_2_btn() {
    MP_btn();
    if (MP3_2.style.display === "none") {
        MP1.style.display = "none";
        MP2_1.style.display = "none";
        MP2_2.style.display = "none";
        MP3_1.style.display = "none";
        MP3_2.style.display = "block";
        MP4.style.display = "none";
        document.getElementById("MP3_2_btn").style.backgroundColor = "#4CAF50";
        document.getElementById("MP3_2_btn").style.color = "white";
    } else {
        MP3_2.style.display = "none";
        document.getElementById("MP3_2_btn").style.backgroundColor = "";
        document.getElementById("MP3_2_btn").style.color = "";
    }
}
function MP4_btn() {
    MP_btn();
    if (MP4.style.display === "none") {
        MP1.style.display = "none";
        MP2_1.style.display = "none";
        MP2_2.style.display = "none";
        MP3_1.style.display = "none";
        MP3_2.style.display = "none";
        MP4.style.display = "block";
        document.getElementById("MP4_btn").style.backgroundColor = "#4CAF50";
        document.getElementById("MP4_btn").style.color = "white";
    } else {
        MP4.style.display = "none";
        document.getElementById("MP4_btn").style.backgroundColor = "";
        document.getElementById("MP4_btn").style.color = "";
    }
}

// This is a test type para in case na mag ctrl + z