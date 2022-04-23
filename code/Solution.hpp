/**
 * @file Solution.hpp
 *
 * @author Vikas Sheoran (f20180847@hyderabad.bits-pilani.ac.in)
 * @author Ritika Garg (f20180256@hyderabad.bits-pilani.ac.in)
 * @author Akshat Goyal (f20180864@hyderabad.bits-pilani.ac.in)
 * @author Vastav Ratra (f20180903@hyderabad.bits-pilani.ac.in)
 *
 * @brief Implements RNA Secondary Structure Prediction Class via Base Pair
 * @version 0.1
 * @bug No known bugs.
 * @date 2022-04-23
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

/* -- Includes -- */

/* libc includes. */
#include <iostream> /* for cin, cout */
#include <string>   /* for string operations */
#include <vector>   /* for 2D array storage */
#include <stack>    /* for printing indices from string */

/*
 *    std namespace declaration throughout
 */
using namespace std;

/**
 * @brief Structure to store indices of bases to pair from given RNA sequence
 * @struct IndexPair
 *
 * @var IndexPair::beginIdx
 * Member 'beginIdx' contains position of the base to pair from the left
 *
 * @var IndexPair::endIdx
 * Member 'endIdx' contains position of the base to pair from the right
 */
struct IndexPair
{
    /*@{*/
    int beginIdx; /**< Contains position of the base to pair from the left */
    int endIdx;   /**< Contains position of the base to pair from the right */
    /*@}*/
};

/**
 * @brief Implements the <b>RNA Secondary Structure Prediction</b> via dynamic programming approach
 * @class Solution
 *
 */
class Solution
{
    vector<vector<int>> dp; /**< 2D table to store memoized results */
    string out;             /**< string to print and visualize pairing */
    string s;               /**< input RNA molecule sequence */
    int n;                  /**< length of RNA molecule sequence */

public:
    /**
     * @brief Construct a new Solution object
     * Empty constuctor necessary
     */
    Solution() {}
    /**
     * @brief Construct a new Solution object
     * And initiliaze our memo table
     * @param rna_sequence
     */
    Solution(string rna_sequence)
    {
        s = rna_sequence;
        n = rna_sequence.length();
        dp.reserve(n);
        for (int i = 0; i < n; i++)
        {
            dp[i] = vector<int>(n, -1); /**< Initalizing 2D table with -1 value implying value not calculated yet */
        }
    }
    /**
     * @brief Perform algorithm and store results in our 2D memo table
     * @return void
     */
    void compute()
    {
        init();
        dp[0][n - 1] = OPT(0, n - 1);
    }
    /**
     * @brief Get the Max Size from all possible RNA secondary foldings
     * @note run after compute ( all caculations are done)
     *
     * @return int
     */
    int getMaxSize()
    {
        return dp[0][n - 1];
    }
    /**
     * @brief Handles printing to stdout all the information post calculations
     *
     */
    void printSet()
    {
        initOut();           /**< Initalizing our output string with '.'s */
        traceback(0, n - 1); /**< Traceback from 2D table to figure relevant pairings */

        printString();  /**< Print the output string */
        printIndices(); /**< Print the indices of all pairs from maximum size RNA secondary folding */
    }

private:
    /**
     * @brief Initialize 2D table with 0 wherever i >= j-4, because of <b> No sharp turn </b> rule
     *
     */
    void init()
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i >= j - 4)
                    dp[i][j] = 0;
            }
        }
    }

    /**
     * @brief Returns the optimal pairings between indices i and j of input RNA sequence
     *
     * @param i starting index
     * @param j ending index
     * @return int maximal size of RNA secondary folding structure set on (i,j)
     */
    int OPT(int i, int j)
    {

        if (i >= j - 4)
            return 0; // No sharp turns
        if (dp[i][j] != -1)
            return dp[i][j]; // Memoization, if already computed return that value

        int exclude = OPT(i, j - 1); // Optimal value when not including current base

        int max_include = 0; // Variable to hold value of optimal pairing when including current base

        // Iterating over all possibilities
        for (int t = i; t < j - 4; t++)
        {
            if (canPair(s[t], s[j])) // Only check if we are pairing valid bases
            {
                int include = 1 + OPT(i, t - 1) + OPT(t + 1, j - 1); // Optimal value when including base upto t
                max_include = max(max_include, include);             // Comparing with maximal value
            }
        }

        return dp[i][j] = max(exclude, max_include); // Store in memo table
    }

    /**
     * @brief Checks whether we are pairing correct bases or not.
     * We can pair only A(adenine) with U(uracil), and C(cytosine) with G(guanine)
     *
     * @param x Base x
     * @param y Base y
     * @return true Yes, we can pair bases x and y
     * @return false No, we can't pair bases x and y
     */
    bool canPair(char x, char y)
    {
        bool ok = (x == 'A' and y == 'U');
        ok |= (x == 'U' && y == 'A');
        ok |= (x == 'C' && y == 'G');
        ok |= (x == 'G' && y == 'C');
        return ok;
    }

    /**
     * @brief Initialize output string with all '.'s.
     * And modify later accounting for all the pairing
     *
     */
    void initOut()
    {
        for (int i = 0; i < n; i++)
        {
            out += ".";
        }
    }

    /**
     * @brief Compute and mark in output string, the optimal pairing by looking at 2D memo table b/w indices i and j
     *
     * @param i start index
     * @param j end index
     */
    void traceback(int i, int j)
    {
        if (i >= j)
            return;                   // Out of bounds
        if (dp[i][j] == dp[i][j - 1]) // Check if we need to include current or skip
        {
            traceback(i, j - 1); // if skip, traceback to sub-table (i,j-1)
            return;
        }

        for (int k = i; k < j; k++) // If can't skip, look at all possibilites
        {
            if (canPair(s[k], s[j])) // Only check if we can pair the bases
            {
                int score = 0;
                if (k - 1 >= 0)
                    score = dp[i][k - 1];
                if (dp[i][j] == 1 + score + dp[k + 1][j - 1]) // If this is the pairing
                {
                    out[k] = '(';            // write to output string startIdx
                    out[j] = ')';            // write to output string endIdx
                    traceback(i, k - 1);     // trackback b/w (i,k-1) now
                    traceback(k + 1, j - 1); // trackback b/w (k+1,j-1) now
                    return;
                }
            }
        }
    }
    /**
     * @brief utility function to print output string showing valid pairing using '.', '(' and ')' to stdout
     *
     */
    void printString()
    {
        cout << "Output: " << out << endl;
    }

    /**
     * @brief utility function to print pairing by showing indices in format <b>(i, j)</b>
     *
     */
    void printIndices()
    {

        vector<IndexPair> indexPairs; /**< vector to stores indices */
        stack<int> index_stack;       /**< helper stack to compute indices from output string */

        // Compute indices from output string, if '(' push into stack. If ')' pop and add to indexPairs
        for (int i = 0; i < out.length(); i++)
        {
            char ch = out[i]; // get current char
            if (ch == '(')
            {
                index_stack.push(i); // push into stack the beginIdx of the base
            }
            else if (ch == ')')
            {
                IndexPair ip; // make a IndexPair struct

                ip.beginIdx = index_stack.top(); // store beginIdx of the base
                ip.endIdx = i;                   // store endIdx of the base
                index_stack.pop();               // pop from stack, since matched already

                indexPairs.push_back(ip); // push into the indexPairs for futher calcuations
            }
        }
        printIndices(indexPairs); // call for printing
    }

    /**
     * @brief utility function to print index pairing from a vector of pairs in format <b> (i, j) </b>
     *
     * @param indexPairs
     */
    void printIndices(vector<IndexPair> const &indexPairs)
    {
        cout << "Indices :" << endl;
        for (IndexPair ip : indexPairs)
        {
            cout << '(' << ip.beginIdx + 1 << ", " << ip.endIdx + 1 << ")" << endl;
        }
    }
};