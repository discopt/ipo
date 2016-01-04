#include "min_norm_2d.h"

#include <cassert>
#include <set>
#include <algorithm>
#include <vector>

#include "vector_2d.h"
#include "rows.h"

#include <iostream> // TODO: Debug

using namespace soplex;

namespace ipo {

  void convergents(const mpz_class& a, const mpz_class& b, std::vector<mpq_class>& result)
  {
    result.clear();
    mpz_class prev = a;
    mpz_class last = b;
    mpz_class matrix11 = 1;
    mpz_class matrix21 = 0;
    mpz_class matrix12 = 0;
    mpz_class matrix22 = 1;
    mpz_class q, r;
    while (last != 0)
    {
      mpz_divmod(q.get_mpz_t(), r.get_mpz_t(), prev.get_mpz_t(), last.get_mpz_t());

      /*
       * n: new matrix, o: old matrix
       * n11 = o11 * q + o12
       * n12 = o11
       * n21 = o21 * q + o22
       * n22 = o21
       */
      std::swap(matrix11, matrix12);
      matrix11 += matrix12 * q;
      std::swap(matrix21, matrix22);
      matrix21 += matrix22 * q;
      result.push_back(mpq_class(matrix11, matrix21));
      prev = last;
      last = r;
    }
  }

  void convergents(const mpq_class& x, std::vector<mpq_class>& result)
  {
    convergents(x.get_num(), x.get_den(), result);
  }

  void hermiteNormalFormTwoRows(const std::vector<mpz_class>& first, const std::vector<mpz_class>& second,
      mpz_class& upperLeft, mpz_class& lowerLeft, mpz_class& lowerRight)
  {
    std::vector<mpz_class> upper, lower;
    assert(first.size() == second.size());
    for (std::size_t i = 0; i < first.size(); ++i)
    {
      if (first[i] != 0)
      {
        upper.push_back(first[i]);
        lower.push_back(second[i]);
      }
    }
    for (std::size_t i = 0; i < first.size(); ++i)
    {
      if (first[i] == 0)
        lower.push_back(second[i]);
    }
    while (upper.size() > 1)
    {
//      std::cerr << "upper:";
//      for (std::size_t i = 0; i < upper.size(); ++i)
//        std::cerr << " " << upper[i];
//      std::cerr << "\nlower:";
//      for (std::size_t i = 0; i < lower.size(); ++i)
//        std::cerr << " " << lower[i];
//      std::cerr << std::endl;

      const std::size_t pivot = upper.size() - 2;
      const mpz_class& x = upper[pivot];
      const mpz_class& y = upper.back();
      mpz_class g, s, t;
      mpz_gcdext(g.get_mpz_t(), s.get_mpz_t(), t.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t());
      mpz_class oldLowerPivot = lower[pivot];
      lower[pivot] = s * oldLowerPivot + t * lower[pivot + 1];
      lower[pivot + 1] = (x * lower[pivot + 1] - y * oldLowerPivot) / g;
      upper[pivot] = g;
      upper.pop_back();
    }
    assert(upper.size() == 1);
    while (lower.size() > 2)
    {
//      std::cerr << "upper:";
//      for (std::size_t i = 0; i < upper.size(); ++i)
//        std::cerr << " " << upper[i];
//      std::cerr << "\nlower:";
//      for (std::size_t i = 0; i < lower.size(); ++i)
//        std::cerr << " " << lower[i];
//      std::cerr << std::endl;

      mpz_gcd(lower[lower.size() - 2].get_mpz_t(), lower[lower.size() - 2].get_mpz_t(), lower.back().get_mpz_t());
      lower.pop_back();
    }
//    std::cerr << "upper:";
//    for (std::size_t i = 0; i < upper.size(); ++i)
//      std::cerr << " " << upper[i];
//    std::cerr << "\nlower:";
//    for (std::size_t i = 0; i < lower.size(); ++i)
//      std::cerr << " " << lower[i];
//    std::cerr << std::endl;

    assert(lower.size() == 2);
    if (lower.back() == 0)
      throw std::runtime_error("HNF is singular!");
    mpz_mod(lower.front().get_mpz_t(), lower.front().get_mpz_t(), lower.back().get_mpz_t());
    upperLeft = upper.front();
    lowerLeft = lower.front();
    lowerRight = lower.back();
  }

  mpz_class approximateLatticeWidthIntegralTriangle(const Integer2d& a, const Integer2d& b, Integer2d& direction)
  {
    mpz_class gamma;
    Integer2d rightMultipliers;
    mpz_gcdext(gamma.get_mpz_t(), rightMultipliers.x.get_mpz_t(), rightMultipliers.y.get_mpz_t(), b.x.get_mpz_t(),
        b.y.get_mpz_t());
    Integer2d leftMultipliers(-b.y / gamma, b.x / gamma);
    mpz_class alpha = leftMultipliers * a;
    if (alpha < 0)
    {
      mpz_neg(alpha.get_mpz_t(), alpha.get_mpz_t());
      mpz_neg(leftMultipliers.x.get_mpz_t(), leftMultipliers.x.get_mpz_t());
      mpz_neg(leftMultipliers.y.get_mpz_t(), leftMultipliers.y.get_mpz_t());
    }
    mpz_class beta = rightMultipliers * a;
    mpz_mod(beta.get_mpz_t(), beta.get_mpz_t(), alpha.get_mpz_t());

//    std::cerr << "HNF of rows " << a << " and " << b << " is\n";
//    std::cerr << alpha << " " << beta << "\n0 " << gamma << "\n";
//    std::cerr << "Via multipliers " << leftMultipliers << " and "
//        << rightMultipliers << std::endl;

/// Check combinations of (alpha,0) and (beta,gamma) according to convergents of a/b.
    std::vector<mpq_class> conv;
    convergents(mpq_class(beta, alpha), conv);
    mpz_class bestMaxNorm = alpha;
    direction = leftMultipliers;
    for (std::size_t i = 0; i < conv.size(); ++i)
    {
      //      std::cerr << "Convergent " << conv[i] << std::endl;
      Integer2d dir(-conv[i].get_num() * alpha + conv[i].get_den() * beta, conv[i].get_den() * gamma);
      dir.normalize();
      mpz_class maxNorm = dir.maximumNorm();
      if (maxNorm < bestMaxNorm)
      {
        direction.x = conv[i].get_den() * rightMultipliers.x - conv[i].get_den() * leftMultipliers.x;
        direction.y = conv[i].get_den() * rightMultipliers.y - conv[i].get_den() * leftMultipliers.y;
        direction.normalize();
        bestMaxNorm = maxNorm;
      }
    }

    /// Compute width along direction.
    mpz_class va = a * direction;
    mpz_class vb = b * direction;
    if (va < vb)
      std::swap(va, vb);
    if (vb > 0)
      return va;
    else if (va < 0)
      return -vb;
    else
      return va - vb;
  }

  void intersectLineWithVerticalLine(const Rational2d& first, const Rational2d& second, const mpz_class& x,
      bool& feasible, bool& onBorder, mpq_class& bestMin, mpq_class& bestMax)
  {
    if (first.x == second.x)
    {
      onBorder = true;
      return;
    }
    if (!(((first.x < second.x) && (x < first.x || x > second.x))
        || ((first.x > second.x) && (x > first.x || x < second.x))))
    {
      mpq_class value = first.y + (second.y - first.y) * (x - first.x) / (second.x - first.x);
      if (!feasible || value < bestMin)
        bestMin = value;
      if (!feasible || value > bestMax)
        bestMax = value;
      feasible = true;
    }
  }

  bool integerProgramOriginThinWidth6Triangle(const Rational2d& left, const Rational2d& right,
      const Integer2d& thinDirection, const Integer2d& objective, Integer2d& solution, mpz_class& solutionObjective)
  {
//    std::cerr << "IP-Triangle" << std::endl;
//    std::cerr << "left*thin = "
//        << mpq_class(left.x * thinDirection.x + left.y * thinDirection.y).get_d()
//        << std::endl;
//    std::cerr << "right*thin = "
//        << mpq_class(right.x * thinDirection.x + right.y * thinDirection.y).get_d()
//        << std::endl;
//    std::cerr << "left*obj = "
//        << mpq_class(left.x * objective.x + left.y * objective.y).get_d()
//        << std::endl;
//    std::cerr << "right*obj = "
//        << mpq_class(right.x * objective.x + right.y * objective.y).get_d()
//        << std::endl;

/// Compute a unimodular transformation whose transposed inverse maps direction thinDirection to (1,0).

    mpz_class gcd, U11, U12, U21, U22;
    mpz_gcdext(gcd.get_mpz_t(), U22.get_mpz_t(), U21.get_mpz_t(), thinDirection.x.get_mpz_t(),
        thinDirection.y.get_mpz_t());
    assert(gcd == 1);
    U21 = -U21;
    U11 = thinDirection.x;
    U12 = thinDirection.y;

    /// Transform left and right vector.
    Rational2d transformedLeft(U11 * left.x + U12 * left.y, U21 * left.x + U22 * left.y);
    Rational2d transformedRight(U11 * right.x + U12 * right.y, U21 * right.x + U22 * right.y);

//    std::cerr << "t_left*(1,0) = " << left.x.get_d() << std::endl;
//    std::cerr << "t_right*(1,0) = " << right.x.get_d() << std::endl;

/// Transform objective via inverse of transpose.
    Integer2d transformedObjective(U22 * objective.x - U21 * objective.y, -U12 * objective.x + U11 * objective.y);

//    std::cerr << "t_left*t_obj = "
//        << mpq_class(
//            transformedLeft.x * transformedObjective.x
//                + transformedLeft.y * transformedObjective.y).get_d()
//        << std::endl;
//    std::cerr << "t_right*t_obj = "
//        << mpq_class(
//            transformedRight.x * transformedObjective.x
//                + transformedRight.y * transformedObjective.y).get_d()
//        << std::endl;
//
//    std::cerr << "Transformed rays are " << transformedLeft << " ~ ("
//        << transformedLeft.x.get_d() << "," << transformedLeft.y.get_d()
//        << ") and " << transformedRight << " ~ (" << transformedRight.x.get_d()
//        << "," << transformedRight.y.get_d() << ")" << std::endl;
//    std::cerr << "Transformed objective is " << transformedObjective
//        << std::endl;

    if (transformedLeft.x > transformedRight.x)
      std::swap(transformedLeft, transformedRight);
    mpq_class minimum = transformedLeft.x < 0 ? transformedLeft.x : mpq_class(0);
    mpq_class maximum = transformedRight.x > 0 ? transformedRight.x : mpq_class(0);
    bool improved = false;
    for (int layer = -6; layer < 6; ++layer)
    {
      if (layer < minimum || layer > maximum)
        continue;

//      std::cerr << "Considering layer tx = " << layer << std::endl;

      bool feasible = false;
      bool onBorder = false;
      mpq_class lowest = 0;
      mpq_class highest = 0;
      intersectLineWithVerticalLine(Rational2d(0, 0), transformedLeft, layer, feasible, onBorder, lowest, highest);
      intersectLineWithVerticalLine(Rational2d(0, 0), transformedRight, layer, feasible, onBorder, lowest, highest);
      intersectLineWithVerticalLine(transformedLeft, transformedRight, layer, feasible, onBorder, lowest, highest);
//      if (onBorder)
//        std::cerr << "Intersection is on border of triangle." << std::endl;
      if (onBorder || lowest == highest)
        continue;
//      std::cerr << "Intersecting segment is [" << lowest << "," << highest
//          << "]." << std::endl;

      mpz_class lowestRounded, highestRounded;
      mpz_cdiv_q(lowestRounded.get_mpz_t(), lowest.get_num_mpz_t(), lowest.get_den_mpz_t());
      if (lowestRounded == lowest)
        lowestRounded++;
      mpz_fdiv_q(highestRounded.get_mpz_t(), highest.get_num_mpz_t(), highest.get_den_mpz_t());
      if (highestRounded == highest)
        highestRounded--;

      if (lowestRounded >= highestRounded)
        continue;

//      std::cerr << "Rounded segment is [" << lowestRounded << ","
//          << highestRounded << "]." << std::endl;

      mpz_class lowestObjective = transformedObjective.x * layer + transformedObjective.y * lowestRounded;
      mpz_class highestObjective = transformedObjective.x * layer + transformedObjective.y * highestRounded;

//      std::cerr << "With objectives " << lowestObjective << " ... "
//          << highestObjective << std::endl;

      mpz_class* yImproved = NULL;
      if (lowestObjective < solutionObjective)
      {
        solution.x = U22 * layer - U12 * lowestRounded;
        solution.y = -U21 * layer + U11 * lowestRounded;
        solutionObjective = lowestObjective;
        improved = true;
      }
      if (highestObjective < solutionObjective)
      {
        solution.x = U22 * layer - U12 * highestRounded;
        solution.y = -U21 * layer + U11 * highestRounded;
        solutionObjective = highestObjective;
        improved = true;
      }
    }
    return improved;
  }

  void manhattanNormShortestLatticeCombination(const std::vector<mpz_class>& u, const std::vector<mpz_class>& v,
      mpz_class& optLambda, mpz_class& optMu, mpz_class& optPi)
  {
    std::size_t n = u.size();
    assert(v.size() == n);
    optLambda = 1;
    optMu = 0;
    optPi = 0;

    mpz_class value;

    /// Find quotients and collect corresponding ray solutions.
    std::set<Integer2d> rays;
    rays.insert(Integer2d(0, 1));
    rays.insert(Integer2d(0, -1));

    /// Add trivial ray covering the (due to lattice basis) possibly changed norm of u.
    rays.insert(Integer2d(1, 0));
    for (std::size_t j = 0; j < n; ++j)
    {
      if (u[j] > 0)
        optPi += u[j];
      else
        optPi -= u[j];
    }

    /// Find other rays.
    mpz_class pi = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
      if (v[i] != 0 && u[i] != 0)
      {
        Integer2d ray = v[i] > 0 ? Integer2d(v[i], -u[i]) : Integer2d(-v[i], u[i]);
        ray.normalize();
        //        std::cerr << "Found ray " << ray << std::endl;
        if (rays.insert(ray).second)
        {
          pi = 0;
          for (std::size_t j = 0; j < n; ++j)
          {
            mpz_class value = ray.x * u[j] + ray.y * v[j];
            mpz_abs(value.get_mpz_t(), value.get_mpz_t());
            pi += value;
          }
          if (pi < optPi)
          {
            optLambda = ray.x;
            optMu = ray.y;
            optPi = pi;
          }
        }
      }
    }

    /// Solve integer programs for each interior of a cone spanned by two rays.
    std::set<Integer2d>::const_iterator leftRay = rays.begin();
    std::set<Integer2d>::const_iterator rightRay = rays.begin();
    Integer2d coneObjective, coneThinDirection, coneMultipliers;
    for (++rightRay; rightRay != rays.end(); ++rightRay)
    {
      //      std::cerr << "\nCone spanned by " << *leftRay << " and " << *rightRay
      //          << std::endl;

      /// Find objective function within that cone.
      coneObjective = Integer2d(0, 0);
      for (std::size_t i = 0; i < n; ++i)
      {
        value = (leftRay->x + rightRay->x) * u[i] + (leftRay->y + rightRay->y) * v[i];
        int sgn = mpz_sgn(value.get_mpz_t());
        coneObjective.x += sgn * u[i];
        coneObjective.y += sgn * v[i];
      }
      //      std::cerr << "Cone objective: " << coneObjective << std::endl;

      /// Goal: Scale both rays such that the endpoint has objective value equal to the best upper bound.
      mpz_class leftScalarDenominator = *leftRay * coneObjective;
      mpz_class rightScalarDenominator = *rightRay * coneObjective;
      Integer2d integerScaledLeftRay(leftRay->x * optPi * rightScalarDenominator,
          leftRay->y * optPi * rightScalarDenominator);
      Integer2d integerScaledRightRay(rightRay->x * optPi * leftScalarDenominator,
          rightRay->y * optPi * leftScalarDenominator);

      /// Find approximate width of that triangle.
      mpz_class width = approximateLatticeWidthIntegralTriangle(integerScaledLeftRay, integerScaledRightRay,
          coneThinDirection);

      mpz_class gcdTest;
      mpz_gcd(gcdTest.get_mpz_t(), coneThinDirection.x.get_mpz_t(), coneThinDirection.y.get_mpz_t());
      assert(gcdTest == 1);

      mpq_class bestUpperBoundWidth(width, leftScalarDenominator * rightScalarDenominator);
      bestUpperBoundWidth.canonicalize();
      //      std::cerr << "Width along thin direction " << coneThinDirection << " is "
      //          << width << "/" << (leftScalarDenominator * rightScalarDenominator)
      //          << " ~ " << bestUpperBoundWidth.get_d() << std::endl;

      mpq_class scaling = mpq_class(5) / width;
      //      mpq_class maxScaling(1, leftScalarDenominator * rightScalarDenominator);
      //      if (scaling > maxScaling)
      //        scaling = maxScaling;

      Rational2d scaledLeftRay(integerScaledLeftRay.x * scaling, integerScaledLeftRay.y * scaling);
      Rational2d scaledRightRay(integerScaledRightRay.x * scaling, integerScaledRightRay.y * scaling);
      //      std::cerr << "Product of scaled left with thin dir is "
      //          << mpq_class(
      //              scaledLeftRay.x * coneThinDirection.x
      //                  + scaledLeftRay.y * coneThinDirection.y).get_d() << std::endl;
      //      std::cerr << "Product of scaled right with thin dir is "
      //          << mpq_class(
      //              scaledRightRay.x * coneThinDirection.x
      //                  + scaledRightRay.y * coneThinDirection.y).get_d()
      //          << std::endl;

      Integer2d multipliers;
      if (integerProgramOriginThinWidth6Triangle(scaledLeftRay, scaledRightRay, coneThinDirection, coneObjective,
          multipliers, optPi))
      {
        optLambda = multipliers.x;
        optMu = multipliers.y;
      }

      //      if (improved)
      //      {
      //        std::cerr << "Found a better solution: " << bestMultipliers
      //            << " with objective " << bestObjective << std::endl;
      //      }

      leftRay = rightRay;
    }

    /// Output best solution.
    //    std::cerr << "Optimal multipliers: " << bestMultipliers << std::endl;
    if (optLambda <= 0)
      throw std::runtime_error("Linear combination has lambda <= 0!");
  }

  void manhattanNormShortestCombination(const std::vector<mpz_class>& u, const std::vector<mpz_class>& v,
      mpq_class& optLambda, mpq_class& optMu, mpz_class& optPi)
  {
    /// Make a quick check whether u and v are a lattice basis already.
    bool uHasUnit = false;
    bool uHasOne = false;
    bool vHasUnit = false;
    bool vHasOne = false;
    std::size_t n = u.size();
    assert(n == v.size());
    mpz_class uNorm = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
      if (!vHasUnit)
      {
        if (v[i] == 1 || v[i] == -1)
        {
          if (u[i] == 0)
            vHasUnit = true;
          vHasOne = true;
        }
      }
      if (!uHasUnit)
      {
        if (u[i] == 1 || u[i] == -1)
        {
          if (v[i] == 0)
            uHasUnit = true;
          uHasOne = true;
        }
      }
      if (u[i] >= 0)
        uNorm += u[i];
      else
        uNorm -= u[i];
    }
    optLambda = 0;
    optMu = 0;
    if ((uHasUnit && vHasOne) || (uHasOne && vHasUnit))
      return manhattanNormShortestLatticeCombination(u, v, optLambda.get_num(), optMu.get_num(), optPi);

    /// If not obviously already a lattice basis, compute HNF.
    mpz_class alpha, beta, gamma;
    hermiteNormalFormTwoRows(v, u, alpha, beta, gamma);
    if (alpha == 1 && gamma == 1)
      return manhattanNormShortestLatticeCombination(u, v, optLambda.get_num(), optMu.get_num(), optPi);

    //    std::cerr << "(Lattice basis reduction with " << upperLeft << ", "
    //        << upperRight << ", " << lowerRight << ")" << std::endl;

    /// It is no lattice basis, so make one!
    std::vector<mpz_class> bu, bv;
    bu.resize(n);
    bv.resize(n);
    mpz_class determinant = alpha * gamma;
    for (std::size_t i = 0; i < n; ++i)
    {
      bu[i] = (alpha * u[i] - beta * v[i]) / determinant;
      bv[i] = v[i] / alpha;
    }
    mpz_class lambda, mu;
    manhattanNormShortestLatticeCombination(bu, bv, lambda, mu, optPi);
    optLambda = mpq_class(lambda) / gamma;
    optMu = (mpq_class(mu) - mpq_class(beta, gamma) * lambda) / alpha;
  }

  bool manhattanNormShortestCombination(std::size_t n, soplex::DSVectorRational& newTarget,
      const soplex::SVectorRational& target, const soplex::SVectorRational& source, soplex::Rational& targetMultiplier,
      soplex::Rational& sourceMultiplier, soplex::Rational& norm)
  {
    newTarget.clear();
    std::vector<mpz_class> u(n);
    std::vector<mpz_class> v(n);
    Rational targetNorm = 0;
    for (int p = target.size() - 1; p >= 0; --p)
    {
      assert(rational2mpzDen(target.value(p)) == 1);
      const Rational& x = target.value(p);
      u[target.index(p)] = rational2mpzNum(x);
      if (x > 0)
        targetNorm += x;
      else
        targetNorm -= x;
    }
    for (int p = source.size() - 1; p >= 0; --p)
    {
      assert(rational2mpzDen(source.value(p)) == 1);
      v[source.index(p)] = rational2mpzNum(source.value(p));
    }

    mpq_class lambda, mu;
    mpz_class pi;
    manhattanNormShortestCombination(u, v, lambda, mu, pi);
    bool result = mpz2rational(pi) < targetNorm;
    if (result)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        mpq_class x = lambda * u[i] + mu * v[i];
        x.canonicalize();
        if (x != 0)
          newTarget.add(i, mpq2rational(x));
      }
      targetMultiplier = mpq2rational(lambda);
      sourceMultiplier = mpq2rational(mu);
      norm = mpq2rational(pi);
    }
    return result;
  }

  bool manhattanNormGreedyCombination(const std::vector<mpz_class>& u, const std::vector<mpz_class>& v,
      mpq_class& lambda, mpq_class& mu, mpq_class& pi)
  {
    std::size_t n = u.size();
    assert(v.size() == n);
    std::vector<bool> processed(n, false);
    std::vector<mpz_class> combination(n, 0);

    /// Init result with u.
    lambda = 1;
    mu = 0;
    pi = 0;
    bool result = false;
    for (std::size_t i = 0; i < n; ++i)
    {
      if (u[i] > 0)
        pi += u[i];
      else
        pi -= u[i];
    }

    /// Go through all coordinates and test how to make them zero in a combination.
    for (std::size_t p = 0; p < n; ++p)
    {
      if (processed[p] || u[p] == 0 || v[p] == 0)
        continue;

      /// Compute combination -u_p * v + v_p * u.
      mpz_class gcd = 0;
      mpz_class norm = 0;
      for (std::size_t i = 0; i < n; ++i)
      {
        combination[i] = v[p] * u[i] - u[p] * v[i];
//        std::cout << "combination#" << i << " = " << combination[i] << ", u#" << i << " = " << u[i] << ", v#" << i
//            << " = " << v[i] << std::endl;
        if (combination[i] == 0)
          processed[i] = true;
        else
        {
          mpz_gcd(gcd.get_mpz_t(), gcd.get_mpz_t(), combination[i].get_mpz_t());
          if (combination[i] > 0)
            norm += combination[i];
          else
            norm -= combination[i];
        }
      }
//      std::cout << "gcd = " << gcd << std::endl;

      norm /= gcd;

      if (norm < pi)
      {
        lambda = mpq_class(v[p], gcd);
        mu = mpq_class(-u[p], gcd);
        pi = norm;
        result = true;
      }
    }
    return result;
  }

  bool manhattanNormGreedyCombination(std::size_t n, soplex::DSVectorRational& newTarget,
      const soplex::SVectorRational& target, const soplex::SVectorRational& source, soplex::Rational& targetMultiplier,
      soplex::Rational& sourceMultiplier, soplex::Rational& norm)
  {
    newTarget.clear();
    std::vector<mpz_class> u(n);
    std::vector<mpz_class> v(n);
    for (int p = target.size() - 1; p >= 0; --p)
    {
      assert(rational2mpzDen(target.value(p)) == 1);
      u[target.index(p)] = rational2mpzNum(target.value(p));
    }
    for (int p = source.size() - 1; p >= 0; --p)
    {
      assert(rational2mpzDen(source.value(p)) == 1);
      v[source.index(p)] = rational2mpzNum(source.value(p));
    }

    mpq_class lambda, mu, pi;
    bool result = manhattanNormGreedyCombination(u, v, lambda, mu, pi);
    if (result)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        mpq_class x = lambda * u[i] + mu * v[i];
        x.canonicalize();
        if (x != 0)
          newTarget.add(i, mpq2rational(x));
      }
      targetMultiplier = mpq2rational(lambda);
      sourceMultiplier = mpq2rational(mu);
      norm = mpq2rational(pi);
    }
    return result;
  }

  void manhattanNormImproveEquations(std::size_t n, soplex::LPRowSetRational& equations)
  {
    scaleRowsIntegral(equations);
    soplex::DSVectorRational newTarget;
    bool improved = true;
    while (improved)
    {
      improved = false;
      for (int i = 0; i < equations.num(); ++i)
      {
        for (int j = 0; j < equations.num(); ++j)
        {
          if (i == j)
            continue;
          soplex::Rational targetMultiplier, sourceMultiplier, newNorm;
          if (manhattanNormGreedyCombination(n, newTarget, equations.rowVector(i), equations.rowVector(j),
              targetMultiplier, sourceMultiplier, newNorm))
          {
            improved = true;
            equations.xtend(i, newTarget.size());
            equations.rowVector_w(i).clear();
            equations.rowVector_w(i) = newTarget;
            equations.lhs_w(i) = targetMultiplier * equations.lhs(i) + sourceMultiplier * equations.lhs(j);
            equations.rhs_w(i) = targetMultiplier * equations.rhs(i) + sourceMultiplier * equations.rhs(j);
          }
        }
      }
    }
    improved = true;
    while (improved)
    {
      improved = false;
      for (int i = 0; i < equations.num(); ++i)
      {
        for (int j = 0; j < equations.num(); ++j)
        {
          if (i == j)
            continue;
          soplex::Rational targetMultiplier, sourceMultiplier, newNorm;
          if (manhattanNormShortestCombination(n, newTarget, equations.rowVector(i), equations.rowVector(j),
              targetMultiplier, sourceMultiplier, newNorm))
          {
            improved = true;
            equations.xtend(i, newTarget.size());
            equations.rowVector_w(i).clear();
            equations.rowVector_w(i) = newTarget;
            equations.lhs_w(i) = targetMultiplier * equations.lhs(i) + sourceMultiplier * equations.lhs(j);
            equations.rhs_w(i) = targetMultiplier * equations.rhs(i) + sourceMultiplier * equations.rhs(j);
          }
        }
      }
    }
  }

  void manhattanNormImproveInequality(std::size_t n, soplex::LPRowRational& inequality,
      const soplex::LPRowSetRational& equations)
  {
    scaleRowIntegral(inequality);
    soplex::DSVectorRational newTarget;
    bool improved = true;
    while (improved)
    {
      improved = false;
      for (int i = 0; i < equations.num(); ++i)
      {
        soplex::Rational targetMultiplier, sourceMultiplier, newNorm;
        if (manhattanNormShortestCombination(n, newTarget, inequality.rowVector(), equations.rowVector(i),
            targetMultiplier, sourceMultiplier, newNorm))
        {
          improved = true;
          inequality.setRowVector(newTarget);
          if (targetMultiplier <= 0)
            throw std::runtime_error("BUG in manhattan norm inequality improvement: Multiplier must be positive!");
          inequality.setLhs(targetMultiplier * inequality.lhs() + sourceMultiplier * equations.lhs(i));
          inequality.setRhs(targetMultiplier * inequality.rhs() + sourceMultiplier * equations.rhs(i));
        }
      }
    }
  }

}
/* namespace ipo */
