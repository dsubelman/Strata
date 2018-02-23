/*
 * Copyright (C) 2018 - present by OpenGamma Inc. and the OpenGamma group of companies
 *
 * Please see distribution for license.
 */
package com.opengamma.strata.market.curve.interpolator;

import java.io.Serializable;

import com.google.common.base.Supplier;
import com.google.common.base.Suppliers;
import com.opengamma.strata.collect.DoubleArrayMath;
import com.opengamma.strata.collect.array.DoubleArray;

/**
 * Extrapolator implementation that is designed for extrapolating a zero rate curve for the near end.
 * <p>
 * The extrapolation is completed by applying a quadratic extrapolant on the discount
 * factor, where the point (0,1) is inserted and
 * the first derivative value is assumed to be continuous at the first x-value.
 */
public class QuadraticLeftZeroRateCurveExtrapolator
    implements CurveExtrapolator, Serializable {

  /**
   * The serialization version id.
   */
  private static final long serialVersionUID = 1L;
  /**
   * The extrapolator name.
   */
  public static final String NAME = "QuadraticLeftZeroRate";
  /**
   * The extrapolator instance.
   */
  public static final CurveExtrapolator INSTANCE = new QuadraticLeftZeroRateCurveExtrapolator();
  /**
   * The epsilon value.
   */
  private static final double EPS = 1e-8;

  /**
   * Restricted constructor.
   */
  private QuadraticLeftZeroRateCurveExtrapolator() {
  }

  // resolve instance
  private Object readResolve() {
    return INSTANCE;
  }

  //-------------------------------------------------------------------------
  @Override
  public String getName() {
    return NAME;
  }

  @Override
  public BoundCurveExtrapolator bind(DoubleArray xValues, DoubleArray yValues,
      BoundCurveInterpolator interpolator) {
    return new Bound(xValues, yValues, interpolator);
  }

  //-------------------------------------------------------------------------
  @Override
  public String toString() {
    return NAME;
  }

  //-------------------------------------------------------------------------
  /**
   * Bound extrapolator.
   */
  static class Bound implements BoundCurveExtrapolator {
    private final int nodeCount;
    private final double firstXValue;
    private final double firstYValue;
    private final double lastXValue;
    private final double eps;
    private final double leftQuadCoef;
    private final double leftLinCoef;
    private final Supplier<DoubleArray> leftSens;

    Bound(DoubleArray xValues, DoubleArray yValues, BoundCurveInterpolator interpolator) {
      this.nodeCount = xValues.size();
      this.firstXValue = xValues.get(0);
      this.firstYValue = Math.exp(-firstXValue * yValues.get(0));
      this.lastXValue = xValues.get(nodeCount - 1);
      double gradient = -yValues.get(0) * firstYValue -
          firstXValue * firstYValue * interpolator.firstDerivative(firstXValue);
      this.eps = EPS * (lastXValue - firstXValue);
      this.leftQuadCoef = gradient / firstXValue - (firstYValue - 1d) / firstXValue / firstXValue;
      this.leftLinCoef = -gradient + 2d * (firstYValue - 1d) / firstXValue;
      this.leftSens = Suppliers.memoize(() -> interpolator.parameterSensitivity(firstXValue + eps));
    }

    //-------------------------------------------------------------------------

    // TODO test handling for small t
    // TODO test cont for positive t -> negative t

    @Override
    public double leftExtrapolate(double xValue) {
      if (firstXValue == 0d) {
        throw new IllegalArgumentException("The trivial point at x = 0 is already included");
      }
      if (Math.abs(xValue) < eps) {
        return -leftLinCoef + (leftQuadCoef - leftLinCoef * leftLinCoef) * xValue +
            xValue * xValue *
                (leftLinCoef * leftLinCoef * leftLinCoef / 3d - 2d * leftQuadCoef * leftLinCoef);
      }
      double df = leftQuadCoef * xValue * xValue + leftLinCoef * xValue + 1d;
      return -Math.log(df) / xValue;
    }

    @Override
    public double leftExtrapolateFirstDerivative(double xValue) {
      if (firstXValue == 0d) {
        throw new IllegalArgumentException("The trivial point at x = 0 is already included");
      }
      if (Math.abs(xValue) < eps) {
        return leftQuadCoef - leftLinCoef * leftLinCoef + xValue *
            (2d * leftLinCoef * leftLinCoef * leftLinCoef / 3d - 4d * leftQuadCoef * leftLinCoef);
      }
      double gradDf = 2d * leftQuadCoef * xValue + leftLinCoef;
      double df = leftQuadCoef * xValue * xValue + leftLinCoef * xValue + 1d;
      return Math.log(df) / Math.pow(xValue, 2) - gradDf / (df * xValue);
    }

    @Override
    public DoubleArray leftExtrapolateParameterSensitivity(double xValue) {
      if (firstXValue == 0d) {
        throw new IllegalArgumentException("The trivial point at x = 0 is already included");
      }
      double[] sensiDf = leftSens.get().toArray();
      double xQuad = xValue * xValue;
      if (Math.abs(xValue) < eps) {
        double factor =
            (2d * leftQuadCoef - leftLinCoef * leftLinCoef - 2d * leftLinCoef / firstXValue) *
                xQuad + (2d * leftLinCoef + 1d / firstXValue) * xValue + 1d;
        DoubleArrayMath.mutateByMultiplication(sensiDf, factor);
        sensiDf[0] += (2d * leftLinCoef / firstXValue / firstXValue +
            2d * leftLinCoef * leftLinCoef / firstXValue -
            4d * leftQuadCoef / firstXValue) * xQuad -
            (1d / firstXValue / firstXValue + 4d * leftLinCoef / firstXValue) * xValue -
            2d / firstXValue;
        return DoubleArray.ofUnsafe(sensiDf);
      }
      double df = leftQuadCoef * xQuad + leftLinCoef * xValue + 1d;
      for (int i = 1; i < nodeCount; i++) {
        double tmp = sensiDf[i] * xValue / eps;
        sensiDf[i] = tmp / firstXValue * xValue - tmp;
      }
      double tmp = (sensiDf[0] - 1d) / eps;
      sensiDf[0] = (tmp / firstXValue - 1d / firstXValue / firstXValue) * xQuad +
          (2d / firstXValue - tmp) * xValue;
      return DoubleArray.ofUnsafe(sensiDf).dividedBy(-xValue * df);
    }

//-------------------------------------------------------------------------
    @Override
    public double rightExtrapolate(double xValue) {
      throw new IllegalArgumentException(
          "QuadraticLeftZeroRateCurveExtrapolator cannot be used for right extrapolation");
    }

    @Override
    public double rightExtrapolateFirstDerivative(double xValue) {
      throw new IllegalArgumentException(
          "QuadraticLeftZeroRateCurveExtrapolator cannot be used for right extrapolation");
    }

    @Override
    public DoubleArray rightExtrapolateParameterSensitivity(double xValue) {
      throw new IllegalArgumentException(
          "QuadraticLeftZeroRateCurveExtrapolator cannot be used for right extrapolation");
    }
  }

}
