/**
 * Copyright (C) 2017 - present by OpenGamma Inc. and the OpenGamma group of companies
 *
 * Please see distribution for license.
 */
package com.opengamma.strata.basics.currency;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NoSuchElementException;
import java.util.Set;

import org.joda.beans.Bean;
import org.joda.beans.BeanBuilder;
import org.joda.beans.BeanDefinition;
import org.joda.beans.ImmutableBean;
import org.joda.beans.JodaBeanUtils;
import org.joda.beans.MetaProperty;
import org.joda.beans.Property;
import org.joda.beans.PropertyDefinition;
import org.joda.beans.impl.direct.DirectFieldsBeanBuilder;
import org.joda.beans.impl.direct.DirectMetaBean;
import org.joda.beans.impl.direct.DirectMetaProperty;
import org.joda.beans.impl.direct.DirectMetaPropertyMap;

import com.google.common.collect.ImmutableMap;

/**
 * An array of currency amounts with the same currency.
 * <p>
 * Each currency amount is labeled, for example, by {@code StandardId}.
 * <p>
 * @param <T>  the type of key for each value
 */
@BeanDefinition(builderScope = "private")
public final class SplitCurrencyAmount<T>
    implements FxConvertible<SplitCurrencyAmount<T>>, ImmutableBean, Serializable {

  @PropertyDefinition(validate = "notNull")
  private final Currency currency;

  @PropertyDefinition(validate = "notNull")
  private final ImmutableMap<T, Double> splitValues;

  //-------------------------------------------------------------------------
  /**
   * Obtains an instance from currency and map.
   * 
   * @param currency  the currency of the values
   * @param splitValues  the split values
   * @return the instance
   */
  public static <T> SplitCurrencyAmount<T> of(Currency currency, Map<T, Double> splitValues) {
    return new SplitCurrencyAmount<T>(currency, splitValues);
  }

  //-------------------------------------------------------------------------
  @Override
  public SplitCurrencyAmount<T> convertedTo(Currency resultCurrency, FxRateProvider rateProvider) {
    Map<T, Double> mutable = new HashMap<>();
    for (Entry<T, Double> entry : splitValues.entrySet()) {
      double converted = rateProvider.convert(entry.getValue(), currency, resultCurrency);
      mutable.put(entry.getKey(), converted);
    }
    return SplitCurrencyAmount.of(resultCurrency, mutable);
  }
  
  //------------------------- AUTOGENERATED START -------------------------
  ///CLOVER:OFF
  /**
   * The meta-bean for {@code SplitCurrencyAmount}.
   * @return the meta-bean, not null
   */
  @SuppressWarnings("rawtypes")
  public static SplitCurrencyAmount.Meta meta() {
    return SplitCurrencyAmount.Meta.INSTANCE;
  }

  /**
   * The meta-bean for {@code SplitCurrencyAmount}.
   * @param <R>  the bean's generic type
   * @param cls  the bean's generic type
   * @return the meta-bean, not null
   */
  @SuppressWarnings("unchecked")
  public static <R> SplitCurrencyAmount.Meta<R> metaSplitCurrencyAmount(Class<R> cls) {
    return SplitCurrencyAmount.Meta.INSTANCE;
  }

  static {
    JodaBeanUtils.registerMetaBean(SplitCurrencyAmount.Meta.INSTANCE);
  }

  /**
   * The serialization version id.
   */
  private static final long serialVersionUID = 1L;

  private SplitCurrencyAmount(
      Currency currency,
      Map<T, Double> splitValues) {
    JodaBeanUtils.notNull(currency, "currency");
    JodaBeanUtils.notNull(splitValues, "splitValues");
    this.currency = currency;
    this.splitValues = ImmutableMap.copyOf(splitValues);
  }

  @SuppressWarnings("unchecked")
  @Override
  public SplitCurrencyAmount.Meta<T> metaBean() {
    return SplitCurrencyAmount.Meta.INSTANCE;
  }

  @Override
  public <R> Property<R> property(String propertyName) {
    return metaBean().<R>metaProperty(propertyName).createProperty(this);
  }

  @Override
  public Set<String> propertyNames() {
    return metaBean().metaPropertyMap().keySet();
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the currency.
   * @return the value of the property, not null
   */
  public Currency getCurrency() {
    return currency;
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the splitValues.
   * @return the value of the property, not null
   */
  public ImmutableMap<T, Double> getSplitValues() {
    return splitValues;
  }

  //-----------------------------------------------------------------------
  @Override
  public boolean equals(Object obj) {
    if (obj == this) {
      return true;
    }
    if (obj != null && obj.getClass() == this.getClass()) {
      SplitCurrencyAmount<?> other = (SplitCurrencyAmount<?>) obj;
      return JodaBeanUtils.equal(currency, other.currency) &&
          JodaBeanUtils.equal(splitValues, other.splitValues);
    }
    return false;
  }

  @Override
  public int hashCode() {
    int hash = getClass().hashCode();
    hash = hash * 31 + JodaBeanUtils.hashCode(currency);
    hash = hash * 31 + JodaBeanUtils.hashCode(splitValues);
    return hash;
  }

  @Override
  public String toString() {
    StringBuilder buf = new StringBuilder(96);
    buf.append("SplitCurrencyAmount{");
    buf.append("currency").append('=').append(currency).append(',').append(' ');
    buf.append("splitValues").append('=').append(JodaBeanUtils.toString(splitValues));
    buf.append('}');
    return buf.toString();
  }

  //-----------------------------------------------------------------------
  /**
   * The meta-bean for {@code SplitCurrencyAmount}.
   * @param <T>  the type
   */
  public static final class Meta<T> extends DirectMetaBean {
    /**
     * The singleton instance of the meta-bean.
     */
    @SuppressWarnings("rawtypes")
    static final Meta INSTANCE = new Meta();

    /**
     * The meta-property for the {@code currency} property.
     */
    private final MetaProperty<Currency> currency = DirectMetaProperty.ofImmutable(
        this, "currency", SplitCurrencyAmount.class, Currency.class);
    /**
     * The meta-property for the {@code splitValues} property.
     */
    @SuppressWarnings({"unchecked", "rawtypes" })
    private final MetaProperty<ImmutableMap<T, Double>> splitValues = DirectMetaProperty.ofImmutable(
        this, "splitValues", SplitCurrencyAmount.class, (Class) ImmutableMap.class);
    /**
     * The meta-properties.
     */
    private final Map<String, MetaProperty<?>> metaPropertyMap$ = new DirectMetaPropertyMap(
        this, null,
        "currency",
        "splitValues");

    /**
     * Restricted constructor.
     */
    private Meta() {
    }

    @Override
    protected MetaProperty<?> metaPropertyGet(String propertyName) {
      switch (propertyName.hashCode()) {
        case 575402001:  // currency
          return currency;
        case 2075486428:  // splitValues
          return splitValues;
      }
      return super.metaPropertyGet(propertyName);
    }

    @Override
    public BeanBuilder<? extends SplitCurrencyAmount<T>> builder() {
      return new SplitCurrencyAmount.Builder<T>();
    }

    @SuppressWarnings({"unchecked", "rawtypes" })
    @Override
    public Class<? extends SplitCurrencyAmount<T>> beanType() {
      return (Class) SplitCurrencyAmount.class;
    }

    @Override
    public Map<String, MetaProperty<?>> metaPropertyMap() {
      return metaPropertyMap$;
    }

    //-----------------------------------------------------------------------
    /**
     * The meta-property for the {@code currency} property.
     * @return the meta-property, not null
     */
    public MetaProperty<Currency> currency() {
      return currency;
    }

    /**
     * The meta-property for the {@code splitValues} property.
     * @return the meta-property, not null
     */
    public MetaProperty<ImmutableMap<T, Double>> splitValues() {
      return splitValues;
    }

    //-----------------------------------------------------------------------
    @Override
    protected Object propertyGet(Bean bean, String propertyName, boolean quiet) {
      switch (propertyName.hashCode()) {
        case 575402001:  // currency
          return ((SplitCurrencyAmount<?>) bean).getCurrency();
        case 2075486428:  // splitValues
          return ((SplitCurrencyAmount<?>) bean).getSplitValues();
      }
      return super.propertyGet(bean, propertyName, quiet);
    }

    @Override
    protected void propertySet(Bean bean, String propertyName, Object newValue, boolean quiet) {
      metaProperty(propertyName);
      if (quiet) {
        return;
      }
      throw new UnsupportedOperationException("Property cannot be written: " + propertyName);
    }

  }

  //-----------------------------------------------------------------------
  /**
   * The bean-builder for {@code SplitCurrencyAmount}.
   * @param <T>  the type
   */
  private static final class Builder<T> extends DirectFieldsBeanBuilder<SplitCurrencyAmount<T>> {

    private Currency currency;
    private Map<T, Double> splitValues = ImmutableMap.of();

    /**
     * Restricted constructor.
     */
    private Builder() {
    }

    //-----------------------------------------------------------------------
    @Override
    public Object get(String propertyName) {
      switch (propertyName.hashCode()) {
        case 575402001:  // currency
          return currency;
        case 2075486428:  // splitValues
          return splitValues;
        default:
          throw new NoSuchElementException("Unknown property: " + propertyName);
      }
    }

    @SuppressWarnings("unchecked")
    @Override
    public Builder<T> set(String propertyName, Object newValue) {
      switch (propertyName.hashCode()) {
        case 575402001:  // currency
          this.currency = (Currency) newValue;
          break;
        case 2075486428:  // splitValues
          this.splitValues = (Map<T, Double>) newValue;
          break;
        default:
          throw new NoSuchElementException("Unknown property: " + propertyName);
      }
      return this;
    }

    @Override
    public Builder<T> set(MetaProperty<?> property, Object value) {
      super.set(property, value);
      return this;
    }

    @Override
    public Builder<T> setString(String propertyName, String value) {
      setString(meta().metaProperty(propertyName), value);
      return this;
    }

    @Override
    public Builder<T> setString(MetaProperty<?> property, String value) {
      super.setString(property, value);
      return this;
    }

    @Override
    public Builder<T> setAll(Map<String, ? extends Object> propertyValueMap) {
      super.setAll(propertyValueMap);
      return this;
    }

    @Override
    public SplitCurrencyAmount<T> build() {
      return new SplitCurrencyAmount<T>(
          currency,
          splitValues);
    }

    //-----------------------------------------------------------------------
    @Override
    public String toString() {
      StringBuilder buf = new StringBuilder(96);
      buf.append("SplitCurrencyAmount.Builder{");
      buf.append("currency").append('=').append(JodaBeanUtils.toString(currency)).append(',').append(' ');
      buf.append("splitValues").append('=').append(JodaBeanUtils.toString(splitValues));
      buf.append('}');
      return buf.toString();
    }

  }

  ///CLOVER:ON
  //-------------------------- AUTOGENERATED END --------------------------
}
