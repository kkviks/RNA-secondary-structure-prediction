/**
 * @file main.cpp
 *
 * @author Vikas Sheoran (f20180847@hyderabad.bits-pilani.ac.in)
 * @author Ritika Garg (f20180256@hyderabad.bits-pilani.ac.in)
 * @author Akshat Goyal (f20180864@hyderabad.bits-pilani.ac.in)
 * @author Vastav Ratra (f20180903@hyderabad.bits-pilani.ac.in)
 *
 * @brief This is the driver code containing main()
 *
 * It takes input from the standard input device stdin,
 * a RNA molecule represented as a string using the symbols A(adenine), U(uracil), C(cytosine), and G(guanine).
 * It printout to standard output,
 * the maximum size of RNA secondary structure set, along with which bases to match( via providing their indices)
 *
 * @version 0.1
 * @date 2022-04-23
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <iostream>
#include <string>

#include "Solution.hpp"
using namespace std;

string test_case = "ACCGGUAGU"; // ans => 2

/**
 * @brief Entry point for the program
 * @return int Execution status of the program
 */
int main()
{
    // Taking input information about RNA molecule through stdin
    cout << "========== < INTPUT > ==========\n";

    cout << "Enter the RNA Sequence: ";
    string RNA_seqeuence;
    cin >> RNA_seqeuence;

    cout << "========== < /INPUT > ==========\n";

    /**
     * @brief Creating instance of our Solution class
     * @param RNA_Sequence string representing the RNA sequence
     * @return Solution class instance
     */
    Solution solution(RNA_seqeuence);

    /**
     * @brief Calling driver method of Solution class instance
     *
     * It initalizes the trivial state of our DP solution and compute other states using top-down approach
     */
    solution.compute();

    // Printing the output to stdout
    cout << "\n========== < OUTPUT > ==========\n";

    /**
     * @brief Extracting the maximum size of RNA secondary structure set from solution instance
     * @return int maxSize contains the maximum size of RNA secondary structure set
     */
    int maxSize = solution.getMaxSize();

    // Printing the maximum size to stdout
    cout << "Max size is = " << maxSize << endl;

    /**
     * @brief Printing the set information to stdout
     */
    solution.printSet();

    cout << "========== < /OUTPUT > ==========\n";
}