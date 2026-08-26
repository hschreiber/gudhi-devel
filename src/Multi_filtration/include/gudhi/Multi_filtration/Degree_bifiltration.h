/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber, David Loiseaux
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Degree_bifiltration.h
 * @author Hannah Schreiber, David Loiseaux
 * @brief Contains the @ref Gudhi::multi_filtration::Degree_bifiltration class.
 */

#ifndef MF_DEGREE_BIFILTRATION_H_
#define MF_DEGREE_BIFILTRATION_H_

#include <cstddef>      //std::size_t
#include <cstring>      //memcpy
#include <stdexcept>
#include <type_traits>  //std::is_arithmetic
#include <vector>

#include <boost/iterator/iterator_facade.hpp>

#include <gudhi/Debug_utils.h>
#include <gudhi/serialization_utils.h>
#include <gudhi/Multi_filtration/multi_filtration_utils.h>

namespace Gudhi::multi_filtration {

/**
 * @class Degree_bifiltration Degree_bifiltration.h gudhi/Multi_filtration/Degree_bifiltration.h
 * @ingroup multi_filtration
 *
 * @brief Class encoding the different generators, i.e., apparition times, of a \f$ k \f$-critical
 * \f$\mathbb R^2\f$-filtration value in a filtration of the particular form:
 * \f$ [(v_0,0), (v_1,1), (v_2,2), ..., (v_{k-1},k-1)] \f$, where all pairs \f$ (v_i,i) \f$
 * represent a generator with two parameters.
 * For example: let \f$ d \f$ be the max degree of a vertex in a complex. Let a vertex be \f$ k \f$-critical if it has
 * degree \f$ d - k + 1 \f$ and an edge if one of its end vertices is \f$ k \f$-critical and the other one
 * \f$ j \f$-critical, \f$ j \geq k \f$. If we filter the complex by degree, the filtration value of a
 * \f$ k \f$-critical cell has \f$ k \f$ generators whose second parameter takes values from 0 to \f$ (k - 1) \f$.
 * The first parameter can then be set for example as the standard radius parameter of a Rips filtration.
 * Note that the set of generators does not have to be minimal (contrary to, e.g., @ref Flat_array_filtration),
 * neither ordered lexicographically.
 * Implements the concept @ref StoragePolicy.
 *
 * @details For more documentation of the public interface, see @ref StoragePolicy.
 *
 * The generators are always internally ordered by the second parameter.
 *
 * `std::numeric_limits<Multi_parameter_filtration_value<Degree_bifiltration>, Co>` will behave such that:
 * - `::has_infinity` returns `true` if and only if `Co` is false,
 * - `::has_quiet_NaN` returns `true`,
 * - `::infinity(int)` returns @ref Degree_bifiltration::inf(int) "",
 * - `::minus_infinity(int)` returns @ref Degree_bifiltration::minus_inf(int) "",
 * - `::max(int)` throws if `Co` is true and otherwise returns a @ref Degree_bifiltration with one generator with
 * first parameter 0 and second parameter `std::numeric_limits<T>::max()`,
 * - `::quiet_NaN(int)` returns @ref Degree_bifiltration::nan(int).
 *
 * @tparam T Arithmetic type of an entry of the second parameter of a filtration value. Has to be **signed** and
 * to implement `std::isnan(T)`, `std::numeric_limits<T>::has_quiet_NaN`, `std::numeric_limits<T>::quiet_NaN()`,
 * `std::numeric_limits<T>::has_infinity`, `std::numeric_limits<T>::infinity()` and `std::numeric_limits<T>::max()`.
 * If `std::numeric_limits<T>::has_infinity` returns `false`, a call to `std::numeric_limits<T>::infinity()`
 * can simply throw. Examples are the native types `double`, `float` and `int`.
 */
template <typename T>
class Degree_bifiltration {
 private:
  template <bool IsConst>
  class Generator_iterator
      : public boost::iterator_facade<Generator_iterator<IsConst>, std::conditional_t<IsConst, const T, T>,
                                      boost::random_access_traversal_tag> {
   private:
    using Base = boost::iterator_facade<Generator_iterator<IsConst>, std::conditional_t<IsConst, const T, T>,
                                        boost::random_access_traversal_tag>;

   public:
    using value_type = std::conditional_t<IsConst, const T, T>;
    using pointer_type = value_type *;
    using reference = typename Base::reference;              // value_type&
    using difference_type = typename Base::difference_type;  // std::ptrdiff_t

    Generator_iterator() : ptr_first_(nullptr), second_(), index_(2) {}

    Generator_iterator(value_type &first, T second) : ptr_first_(&first), second_(second), index_(0) {}

    // Allow mutable -> const iterator conversion
    template <bool OtherConst, typename = std::enable_if_t<IsConst && !OtherConst>>
    Generator_iterator(const Generator_iterator<OtherConst> &other)
        : ptr_first_(other.ptr_first_), second_(other.second_), index_(other.index_) {}

    // necessary for boost::iterator_range to be able to use operator[].
    // Seg fails otherwise for some reasons
    reference operator[](difference_type n) const {
      GUDHI_CHECK(n == -1 || n == 0 || n == 1, std::out_of_range("Out of range index for Generator."));
      if (n < 0) return second_;
      return (index_ == 0) ? *ptr_first_ : second_;
    }

   private:
    friend class boost::iterator_core_access;
    template <bool>
    friend class Generator_iterator;

    pointer_type ptr_first_;
    mutable T second_;  // mutable because dereference() has to be const
    difference_type index_;

    reference dereference() const {
      GUDHI_CHECK(index_ == 0 || index_ == 1, std::out_of_range("Out of range index for Generator."));
      return (index_ == 0) ? *ptr_first_ : second_;
    }

    bool equal(const Generator_iterator &other) const {
      if (index_ == 2) return other.index_ == 2;
      return ptr_first_ == other.ptr_first_ && second_ == other.second_ && index_ == other.index_;
    }

    void increment() { ++index_; }

    void decrement() { --index_; }

    void advance(difference_type n) { index_ += n; }

    difference_type distance_to(const Generator_iterator &other) const { return other.index_ - index_; }
  };

 public:
  using value_type = T;                                       /**< `T` */
  using Underlying_container = std::vector<T>;                /**< std::vector<T> */
  using size_type = typename Underlying_container::size_type; /**< std::size_t */
  using reference = value_type &;                             /**< @ref value_type & */
  using const_reference = value_type;                         /**< @ref value_type */
  using iterator = Generator_iterator<false>;                 /**< LegacyRandomAccessIterator Iterator type. */
  using const_iterator = Generator_iterator<true>;            /**< LegacyRandomAccessIterator Const iterator type. */
  template <typename U>
  using As_type = Degree_bifiltration<U>; /**< Degree_bifiltration<U> */

  template <typename U = T>
  constexpr static const U T_inf = detail::MF_T_inf<U>;

  template <typename U = T>
  constexpr static const U T_m_inf = detail::MF_T_m_inf<U>;

  constexpr static const bool has_an_implicit_axis = true;
  constexpr static const bool has_lexicographical_storage = false;
  constexpr static const bool has_minimal_set_representation = false;
  template <bool Co>
  constexpr static const bool has_infinity = !Co;
  template <bool Co>
  constexpr static const bool has_minus_infinity = true;
  constexpr static const bool has_quiet_NaN = true;

  Degree_bifiltration() = default;

  Degree_bifiltration([[maybe_unused]] size_type numberOfParameters, T value) : generators_(1, value) {}

  template <class Iterator, class = std::enable_if_t<!std::is_arithmetic_v<Iterator>>>
  Degree_bifiltration(Iterator itBegin, [[maybe_unused]] Iterator itEnd) : generators_(1, *itBegin) {
    GUDHI_CHECK(*(std::next(itBegin)) == 0, std::invalid_argument("Second value of the range has to be 0"));
  }

  template <class Iterator, class = std::enable_if_t<!std::is_arithmetic_v<Iterator>>>
  Degree_bifiltration(Iterator itBegin, Iterator itEnd, [[maybe_unused]] size_type numberOfParameters) {
    size_type numGen = std::distance(itBegin, itEnd) / 2;
    generators_.resize(numGen);
    Iterator it = itBegin;
    for (size_type i = 0; i < numGen; ++i) {
      generators_[i] = *it;
      ++it;
      GUDHI_CHECK(
          static_cast<size_type>(*it) == i,
          std::invalid_argument(
              "Every second value of the range has to correspond to a contiguous sequence of integers starting at 0."));
      ++it;
    }
  }

  Degree_bifiltration(const Underlying_container &generators, [[maybe_unused]] size_type numberOfParameters)
      : generators_(generators) {}

  Degree_bifiltration(Underlying_container &&generators, [[maybe_unused]] size_type numberOfParameters)
      : generators_(std::move(generators)) {}

  Degree_bifiltration(const Degree_bifiltration &other) = default;

  Degree_bifiltration(Degree_bifiltration &&other) noexcept = default;

  ~Degree_bifiltration() = default;

  // cannot use = default as it triggers "dummy_g_ may be used uninitialized" compiler warning for nothing
  Degree_bifiltration &operator=(const Degree_bifiltration &other) {
    generators_ = other.generators_;
    return *this;
  }

  // cannot use = default as it triggers "dummy_g_ may be used uninitialized" compiler warning for nothing
  Degree_bifiltration &operator=(Degree_bifiltration &&other) noexcept {
    generators_ = std::move(other.generators_);
    return *this;
  }

  friend void swap(Degree_bifiltration &f1, Degree_bifiltration &f2) noexcept { f1.generators_.swap(f2.generators_); }

  Underlying_container &get_underlying_container() { return generators_; }

  const Underlying_container &get_underlying_container() const { return generators_; }

  reference operator()(size_type g, size_type p) {
    GUDHI_CHECK(g < generators_.size(), std::out_of_range("Out of bound generator index."));
    GUDHI_CHECK(p < 2, std::out_of_range("Out of bound parameter index."));
    if (p == 1) {
      dummy_g_ = g;
      return dummy_g_;
    }
    return generators_[g];
  }

  const_reference operator()(size_type g, size_type p) const {
    GUDHI_CHECK(g < generators_.size(), std::out_of_range("Out of bound generator index."));
    GUDHI_CHECK(p < 2, std::out_of_range("Out of bound parameter index."));
    if (p == 1) return g;
    return generators_[g];
  }

  iterator begin(size_type g) { return iterator(generators_[g], g); }

  iterator end(size_type g) { return iterator(); }

  const_iterator begin(size_type g) const { return const_iterator(generators_[g], g); }

  const_iterator end(size_type g) const { return const_iterator(); }

  static constexpr size_type implicit_axis() noexcept { return 1; }

  static int get_generator_of_implicit_value(value_type val) {
    if (detail::_is_nan(val) || val == T_m_inf<>) return null_value<int>();
    if (val == T_inf<>) return T_inf<int>;
    auto g = static_cast<int>(val);
    if (g < 0) return null_value<int>();
    return g;
  }

  template <typename U>
  static constexpr U null_value() noexcept {
    return static_cast<U>(-1);  // equal to T_inf for unsigned integers, but we do not allow them anyway
  }

  void reserve(size_type number_of_generators) { generators_.reserve(number_of_generators); }

  static constexpr size_type num_parameters() { return 2; }

  size_type num_generators() const { return generators_.size(); }

  size_type num_entries() const { return generators_.size() * 2; }

  /**
   * @throw std::logic_error If `Co` is true.
   */
  template <bool Co>
  static Degree_bifiltration inf(size_type numberOfParameters = 2) {
    if constexpr (Co) {
      throw std::logic_error("No biggest value possible for Co-filtrations yet.");
    } else {
      return Degree_bifiltration(numberOfParameters, T_inf<>);
    }
  }

  template <bool Co>
  static Degree_bifiltration minus_inf(size_type numberOfParameters = 2) {
    return Degree_bifiltration(numberOfParameters, T_m_inf<>);
  }

  static Degree_bifiltration nan(size_type numberOfParameters = 2) { return Degree_bifiltration(); }

  template <bool Co>
  [[nodiscard]] constexpr bool is_plus_inf() const {
    if constexpr (Co) {
      return false;
    } else {
      if (generators_.empty()) return false;
      for (const T &v : generators_) {
        if (v != T_inf<>) return false;
      }
      return true;
    }
  }

  template <bool Co>
  [[nodiscard]] bool is_minus_inf() const {
    if constexpr (Co) {
      return generators_.size() == 1 && generators_[0] == T_m_inf<>;
    } else {
      return !generators_.empty() && generators_[0] == T_m_inf<>;
    }
  }

  [[nodiscard]] bool is_nan() const { return generators_.empty(); }

  template <bool Co>
  [[nodiscard]] bool is_finite() const {
    if constexpr (Co) {
      return !generators_.empty() && (generators_.size() != 1 || generators_[0] != T_m_inf<>);
    } else {
      if (generators_.empty() || generators_[0] == T_m_inf<>) return false;
      for (const T &v : generators_) {
        if (v != T_inf<>) return true;
      }
      return false;
    }
  }

  void set_num_generators(size_type g, value_type defaultValue) { generators_.resize(g, defaultValue); }

  template <bool Co, class Iterator>
  bool emplace(Iterator startIt, Iterator endIt, value_type nullValue) {
    GUDHI_CHECK(std::distance(startIt, endIt) == 2,
                std::invalid_argument("Wrong range size. Should correspond to the number of parameters."));

    const T val = *startIt;
    ++startIt;

    GUDHI_CHECK(*startIt >= 0, std::invalid_argument("Second parameter has to be a positive index."));

    const size_type index = *startIt;

    if (detail::_is_nan(val)) return false;

    if (index < generators_.size()) {
      if (_dominates<Co>(val, generators_[index])) return false;
      generators_[index] = val;
      return true;
    }

    generators_.resize(index + 1, nullValue);
    generators_[index] = val;
    return true;
  }

  void sort() {}

  /**
   * @brief Outstream operator.
   */
  friend std::ostream &operator<<(std::ostream &stream, const Degree_bifiltration &f) {
    const size_type num_gen = f.num_generators();
    const size_type num_param = f.num_parameters();

    stream << "( k = " << num_gen << " ) ( p = " << num_param << " ) [ ";
    for (size_type i = 0; i < f.generators_.size(); ++i) {
      stream << f.generators_[i] << " ";
    }
    stream << "]";

    return stream;
  }

  friend char *serialize_value_to_char_buffer(const Degree_bifiltration &value, char *start) {
    return serialize_value_to_char_buffer(value.generators_, start);
  }

  friend const char *deserialize_value_from_char_buffer(Degree_bifiltration &value, const char *start) {
    return deserialize_value_from_char_buffer(value.generators_, start);
  }

  friend std::size_t get_serialization_size_of(const Degree_bifiltration &value) {
    return get_serialization_size_of(value.generators_);
  }

 private:
  Underlying_container generators_; /**< Container of the filtration value elements. */
  T dummy_g_;                       /**< Dummy value for the first parameter, such that a reference can be returned. */

  template <bool Co>
  static bool _dominates(T a, T b) {
    if constexpr (Co) {
      return a <= b;
    } else {
      return a >= b;
    }
  }
};

}  // namespace Gudhi::multi_filtration

#endif  // MF_DEGREE_BIFILTRATION_H_
