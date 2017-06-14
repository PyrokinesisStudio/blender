/*
* Copyright 2014 Google Inc. All rights reserved.
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef MATHFU_VECTOR_H_
#define MATHFU_VECTOR_H_

#include "mathfu/utilities.h"

#include <math.h>

/// @file mathfu/vector.h Vector
/// @brief Vector class and functions.
/// @addtogroup mathfu_vector

// Disable spurious warnings generated by MATHFU_VECTOR_OPERATION().
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4127)  // conditional expression is constant
#if _MSC_VER >= 1900             // MSVC 2015
#pragma warning(disable : 4456)  // allow shadowing in unrolled loops
#endif                           // _MSC_VER >= 1900
#elif defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Warray-bounds"
#endif

/// @cond MATHFU_INTERNAL
#define MATHFU_VECTOR_OPERATION(OP) MATHFU_UNROLLED_LOOP(i, d, OP)
/// @endcond

/// @cond MATHFU_INTERNAL
#define MATHFU_VECTOR_OPERATOR(OP)           \
  {                                          \
    Vector<T, d> result;                     \
    MATHFU_VECTOR_OPERATION(result[i] = OP); \
    return result;                           \
  }
/// @endcond

/// @cond MATHFU_INTERNAL
#define MATHFU_VECTOR_SELF_OPERATOR(OP) \
  {                                     \
    MATHFU_VECTOR_OPERATION(OP);        \
    return *this;                       \
  }
/// @endcond

namespace mathfu {

template <class T, int d>
class Vector;

/// @cond MATHFU_INTERNAL
template <class T, int d>
static inline T DotProductHelper(const Vector<T, d>& v1,
                                 const Vector<T, d>& v2);
template <class T>
static inline T DotProductHelper(const Vector<T, 2>& v1,
                                 const Vector<T, 2>& v2);
template <class T>
static inline T DotProductHelper(const Vector<T, 3>& v1,
                                 const Vector<T, 3>& v2);
template <class T>
static inline T DotProductHelper(const Vector<T, 4>& v1,
                                 const Vector<T, 4>& v2);

template <typename T, int d, typename CompatibleT>
static inline Vector<T, d> FromTypeHelper(const CompatibleT& compatible);

template <typename T, int d, typename CompatibleT>
static inline CompatibleT ToTypeHelper(const Vector<T, d>& v);
/// @endcond

/// @addtogroup mathfu_vector
/// @{

/// @class VectorPacked "mathfu/vector.h"
/// @brief Packed N-dimensional vector.
///
/// Some Vector classes are padded so that it's possible to use the data
/// structures with SIMD instructions.  This structure can be used in
/// conjunction with unpacked Vector classes to pack data
/// into flat arrays suitable for sending to a GPU (e.g vertex buffers).
///
/// <p>
/// For example, to pack (store) an unpacked to packed vector:<br>
/// <blockquote><code><pre>
/// VectorPacked<float, 3> packed;
/// Vector<float, 3> vector(3, 2, 1);
/// vector.Pack(&packed);
/// </pre></code></blockquote>
/// or<br>
/// <blockquote><code><pre>
/// Vector<float, 3> vector(3, 2, 1);
/// VectorPacked<float, 3> packed = vector;
/// </pre></code></blockquote>
/// </p>
///
/// <p>
/// To initialize a vector from a packed vector:<br>
/// <blockquote><code><pre>
/// VectorPacked<float, 3> packed = { 3, 2, 1 };
/// Vector<float, 3> vector(packed);
/// </pre></code></blockquote>
///
/// @tparam T type of VectorPacked elements.
/// @tparam d dimensions (number of elements) in the VectorPacked structure.
template <class T, int d>
struct VectorPacked {
  /// Create an uninitialized VectorPacked.
  VectorPacked() {}

  /// Create a VectorPacked from a Vector.
  ///
  /// Both VectorPacked and Vector must have the same number of dimensions.
  /// @param vector Vector to create the VectorPacked from.
  explicit VectorPacked(const Vector<T, d>& vector) { vector.Pack(this); }

  /// Copy a Vector to a VectorPacked.
  ///
  /// Both VectorPacked and Vector must have the same number of dimensions.
  /// @param vector Vector to copy to the VectorPacked.
  /// @returns A reference to this VectorPacked.
  VectorPacked& operator=(const Vector<T, d>& vector) {
    vector.Pack(this);
    return *this;
  }

  /// Elements of the packed vector one per dimension.
  T data[d];
};
/// @}

/// @addtogroup mathfu_vector
/// @{
/// @class Vector "mathfu/vector.h"
/// @brief Vector of d elements with type T
///
/// Vector stores <b>d</b> elements of type <b>T</b> and provides a set
/// functions to perform operations on the set of elements.
///
/// @tparam T type of Vector elements.
/// @tparam d dimensions (number of elements) in the Vector structure.
template <class T, int d>
class Vector {
 public:
  /// @brief Element type to enable reference by other classes.
  typedef T Scalar;

  /// @brief Create an uninitialized Vector.
  inline Vector() {}

  /// @brief Create a vector from another vector copying each element.
  ///
  /// @param v Vector that the data will be copied from.
  inline Vector(const Vector<T, d>& v) {
    MATHFU_VECTOR_OPERATION(data_[i] = v.data_[i]);
  }

  /// @brief Create a vector from another vector of a different type.
  ///
  /// This copies each element of a Vector which makes it possible to between
  /// vectors of different types, for example
  /// <code>float/double/int</code> vectors.
  /// @param v Vector that the data will be copied from.
  /// @tparam U type of Vector elements to copy.
  template <typename U>
  explicit inline Vector(const Vector<U, d>& v) {
    MATHFU_VECTOR_OPERATION(data_[i] = static_cast<T>(v[i]));
  }

  /// @brief Create a vector from a single float.
  ///
  /// Each elements is set to be equal to the value given.
  /// @param s Scalar value that the vector will be initialized to.
  explicit inline Vector(const T& s) { MATHFU_VECTOR_OPERATION(data_[i] = s); }

  /// @brief Create a vector form the first d elements of an array.
  ///
  /// @param a Array of values that the vector will be iniitlized to.
  explicit inline Vector(const T* a) {
    MATHFU_VECTOR_OPERATION(data_[i] = a[i]);
  }

  /// @brief Create a vector from two values.
  ///
  /// @note This method only works when the vector is of size two.
  ///
  /// @param s1 Scalar value for the first element of the vector.
  /// @param s2 Scalar value for the second element of the vector.
  inline Vector(const T& s1, const T& s2) {
    MATHFU_STATIC_ASSERT(d == 2);
    data_[0] = s1;
    data_[1] = s2;
  }

  /// @brief Create a vector from three values.
  ///
  /// @note This method only works when the vector is of size three.
  ///
  /// @param s1 Scalar value for the first element of the vector.
  /// @param s2 Scalar value for the second element of the vector.
  /// @param s3 Scalar value for the third element of the vector.
  inline Vector(const T& s1, const T& s2, const T& s3) {
    MATHFU_STATIC_ASSERT(d == 3);
    data_[0] = s1;
    data_[1] = s2;
    data_[2] = s3;
  }

  /// @brief Create a vector from a 2 component vector and a third value.
  ///
  /// @note This method only works when the vector is of size three.
  ///
  /// @param v12 Vector containing the first 2 values.
  /// @param s3 Scalar value for the third element of the vector.
  inline Vector(const Vector<T, 2>& v12, const T& s3) {
    MATHFU_STATIC_ASSERT(d == 3);
    data_[0] = v12[0];
    data_[1] = v12[1];
    data_[2] = s3;
  }

  /// @brief Create a vector from four values.
  ///
  /// @note This method only works when the vector is of size four.
  ///
  /// @param s1 Scalar value for the first element of the vector.
  /// @param s2 Scalar value for the second element of the vector.
  /// @param s3 Scalar value for the third element of the vector.
  /// @param s4 Scalar value for the forth element of the vector.
  inline Vector(const T& s1, const T& s2, const T& s3, const T& s4) {
    MATHFU_STATIC_ASSERT(d == 4);
    data_[0] = s1;
    data_[1] = s2;
    data_[2] = s3;
    data_[3] = s4;
  }

  /// @brief Create a 4-dimensional vector from a Vector<T, 3>.
  ///
  /// The last element is initialized to the specified value.
  /// @note This method only works with 4 element vectors.
  ///
  /// @param vector3 Vector used to initialize the first 3 elements.
  /// @param value Value used to set the last element of the vector.
  inline Vector(const Vector<T, 3>& vector3, const T& value) {
    MATHFU_STATIC_ASSERT(d == 4);
    data_[0] = vector3[0];
    data_[1] = vector3[1];
    data_[2] = vector3[2];
    data_[3] = value;
  }

  /// @brief Create a vector from two 2 component vectors.
  ///
  /// @note This method only works when the vector is of size four.
  ///
  /// @param v12 Vector containing the first 2 values.
  /// @param v34 Vector containing the last 2 values.
  inline Vector(const Vector<T, 2>& v12, const Vector<T, 2>& v34) {
    MATHFU_STATIC_ASSERT(d == 4);
    data_[0] = v12[0];
    data_[1] = v12[1];
    data_[2] = v34[0];
    data_[3] = v34[1];
  }

  /// @brief Create a vector from packed vector (VectorPacked).
  ///
  /// @param vector Packed vector used to initialize an unpacked.
  explicit inline Vector(const VectorPacked<T, d>& vector) {
    MATHFU_VECTOR_OPERATION(data_[i] = vector.data[i]);
  }

  /// @brief Access an element of the vector.
  ///
  /// @param i Index of the element to access.
  /// @return A reference to the accessed data that can be modified by the
  /// caller.
  inline T& operator()(const int i) { return data_[i]; }

  /// @brief Access an element of the vector.
  ///
  /// @param i Index of the element to access.
  /// @return A reference to the accessed data.
  inline const T& operator()(const int i) const { return data_[i]; }

  /// @brief Access an element of the vector.
  ///
  /// @param i Index of the element to access.
  /// @return A reference to the accessed data that can be modified by the
  /// caller.
  inline T& operator[](const int i) { return data_[i]; }

  /// @brief Access an element of the vector.
  ///
  /// @param i Index of the element to access.
  /// @return A const reference to the accessed.
  inline const T& operator[](const int i) const { return data_[i]; }

  /// @brief GLSL style 3 element accessor.
  ///
  /// This only works with vectors that contain more than 3 elements.
  /// @returns A 3-dimensional Vector containing the first 3 elements of
  // this Vector.
  inline Vector<T, 3> xyz() {
    MATHFU_STATIC_ASSERT(d > 3);
    return Vector<T, 3>(data_[0], data_[1], data_[2]);
  }

  /// @brief GLSL style 3 element accessor.
  ///
  /// This only works with vectors that contain more than 3 elements.
  /// @returns A 3-dimensional Vector containing the first 3 elements of
  // this Vector.
  inline const Vector<T, 3> xyz() const {
    MATHFU_STATIC_ASSERT(d > 3);
    return Vector<T, 3>(data_[0], data_[1], data_[2]);
  }

  /// @brief GLSL style 2 element accessor.
  ///
  /// This only works with vectors that contain more than 2 elements.
  /// @returns A 2-dimensional Vector with the first 2 elements of this Vector.
  inline Vector<T, 2> xy() {
    MATHFU_STATIC_ASSERT(d > 2);
    return Vector<T, 2>(data_[0], data_[1]);
  }

  /// @brief GLSL style 2 element accessor.
  ///
  /// This only works with vectors that contain more than 2 elements.
  /// @returns A 2-dimensional Vector with the first 2 elements of this Vector.
  inline const Vector<T, 2> xy() const {
    MATHFU_STATIC_ASSERT(d > 2);
    return Vector<T, 2>(data_[0], data_[1]);
  }

  /// @brief GLSL style 2 element accessor.
  ///
  /// This only works with vectors that contain 4 elements.
  /// @returns A 2-dimensional Vector with the last 2 elements of this Vector.
  inline Vector<T, 2> zw() {
    MATHFU_STATIC_ASSERT(d == 4);
    return Vector<T, 2>(data_[2], data_[3]);
  }

  /// @brief GLSL style 2 element accessor.
  ///
  /// This only works with vectors that contain 4 elements.
  /// @returns A 2-dimensional Vector with the last 2 elements of this Vector.
  inline const Vector<T, 2> zw() const {
    MATHFU_STATIC_ASSERT(d == 4);
    return Vector<T, 2>(data_[2], data_[3]);
  }

  /// @brief Pack a Vector to a packed "d" element vector structure.
  ///
  /// @param vector Packed "d" element vector to write to.
  inline void Pack(VectorPacked<T, d>* const vector) const {
    MATHFU_VECTOR_OPERATION(vector->data[i] = data_[i]);
  }

  /// @brief Pack a Vector to a array "d" element.
  ///
  /// @param vector array "d" element to write to.
  inline void Pack(T *a) const {
    MATHFU_VECTOR_OPERATION(a[i] = data_[i]);
  }

  /// @brief Return the array of this vector.
  ///
  /// @return The array of this vector.
  inline const T* const Data() const {
    return data_;
  }

  /// @brief Calculate the squared length of this vector.
  ///
  /// @return The length of this vector squared.
  inline T LengthSquared() const { return LengthSquaredHelper(*this); }

  /// @brief Calculate the length of this vector.
  ///
  /// @return The length of this vector.
  inline T Length() const { return LengthHelper(*this); }

  /// @brief Normalize this vector in-place.
  ///
  /// @return The length of this vector.
  inline T Normalize() { return NormalizeHelper(*this); }

  /// @brief Calculate the normalized version of this vector.
  ///
  /// @return The normalized vector.
  inline Vector<T, d> Normalized() const { return NormalizedHelper(*this); }

  static inline bool FuzzyZero(const Vector<T, d>& v) {
    return FuzzyZeroHelper(v);
  }

  /// @brief Load from any type that is some formulation of a length d array of
  ///        type T.
  ///
  /// Essentially this is just a type cast and a load, but it happens safely
  /// so that we avoid aliasing bugs.
  ///
  /// @return `compatible` cast to `Vector<T,d>` and dereferenced.
  template <typename CompatibleT>
  static inline Vector<T, d> FromType(const CompatibleT& compatible) {
    return FromTypeHelper<T, d, CompatibleT>(compatible);
  }

  /// @brief Load into any type that is some formulation of a length d array of
  ///        type T.
  ///
  /// Essentially this is just a type cast and a load, but it happens safely
  /// so that we avoid aliasing bugs.
  ///
  /// @return `v` cast to `CompatibleT` and dereferenced.
  template <typename CompatibleT>
  static inline CompatibleT ToType(const Vector<T, d>& v) {
    return ToTypeHelper<T, d, CompatibleT>(v);
  }

  /// @brief Calculate the dot product of two vectors.
  ///
  /// @param v1 First vector.
  /// @param v2 Second vector.
  /// @return The dot product of v1 and v2.
  static inline T DotProduct(const Vector<T, d>& v1, const Vector<T, d>& v2) {
    return DotProductHelper(v1, v2);
  }

  /// @brief Calculate the hadamard or componentwise product of two vectors.
  ///
  /// @param v1 First vector.
  /// @param v2 Second vector.
  /// @return The hadamard product of v1 and v2.
  static inline Vector<T, d> HadamardProduct(const Vector<T, d>& v1,
                                             const Vector<T, d>& v2) {
    return HadamardProductHelper(v1, v2);
  }

  /// @brief Calculate the cross product of two vectors.
  ///
  /// Note that this function is only defined for 3-dimensional Vectors.
  /// @param v1 First vector.
  /// @param v2 Second vector.
  /// @return The cross product of v1 and v2.
  static inline Vector<T, 3> CrossProduct(const Vector<T, 3>& v1,
                                          const Vector<T, 3>& v2) {
    return CrossProductHelper(v1, v2);
  }

  /// @brief Linearly interpolate two vectors.
  ///
  /// @param v1 First vector.
  /// @param v2 Second vector.
  /// @param percent Percentage from v1 to v2 in range 0.0...1.0.
  /// @return The hadamard product of v1 and v2.
  static inline Vector<T, d> Lerp(const Vector<T, d>& v1,
                                  const Vector<T, d>& v2, const T percent) {
    return LerpHelper(v1, v2, percent);
  }

  /// @brief Generates a random vector.
  ///
  /// The range of each component is bounded by min and max.
  /// @param min Minimum value of the vector.
  /// @param max Maximum value of the vector.
  static inline Vector<T, d> RandomInRange(const Vector<T, d>& min,
                                           const Vector<T, d>& max) {
    return RandomInRangeHelper(min, max);
  }

  /// @brief Compare each component and returns max values.
  ///
  /// @param v1 First vector.
  /// @param v2 Second vector.
  /// @return Max value of v1 and v2.
  static inline Vector<T, d> Max(const Vector<T, d>& v1,
                                 const Vector<T, d>& v2) {
    return MaxHelper(v1, v2);
  }

  /// @brief Compare each component and returns min values.
  ///
  /// @param v1 First vector.
  /// @param v2 Second vector.
  /// @return Min value of v1 and v2.
  static inline Vector<T, d> Min(const Vector<T, d>& v1,
                                 const Vector<T, d>& v2) {
    return MinHelper(v1, v2);
  }

  MATHFU_DEFINE_CLASS_SIMD_AWARE_NEW_DELETE

  /// Elements of the vector.
  T data_[d];
};
/// @}

/// @addtogroup mathfu_vector
/// @{

/// @brief Compare 2 Vectors of the same size for equality.
///
/// @note: The likelyhood of two float values being the same is very small.
/// Instead consider comparing the difference between two float vectors using
/// LengthSquared() with an epsilon value.
/// For example, v1.LengthSquared(v2) < epsilon.
///
/// @return true if the 2 vectors contains the same value, false otherwise.
template <class T, int d>
inline bool operator==(const Vector<T, d>& lhs, const Vector<T, d>& rhs) {
  for (int i = 0; i < d; ++i) {
    if (lhs[i] != rhs[i]) return false;
  }
  return true;
}

/// @brief Compare 2 Vectors of the same size for inequality.
///
/// @return true if the elements of two vectors differ, false otherwise.
template <class T, int d>
inline bool operator!=(const Vector<T, d>& lhs, const Vector<T, d>& rhs) {
  return !(lhs == rhs);
}

/// @brief Negate all elements of the Vector.
///
/// @return A new Vector containing the result.
template <class T, int d>
inline Vector<T, d> operator-(const Vector<T, d>& v) {
  MATHFU_VECTOR_OPERATOR(-v.data_[i]);
}

/// @brief Multiply a Vector by a scalar.
///
/// Multiplies each component of the specified Vector with a scalar.
///
/// @param s scalar to multiply.
/// @param v Vector to multiply.
/// @return Vector containing the result.
/// @related Vector
template <class T, int d>
inline Vector<T, d> operator*(const T& s, const Vector<T, d>& v) {
  MATHFU_VECTOR_OPERATOR(v.data_[i] * s);
}

/// @brief Divide a Vector by a scalar.
///
/// Divides each component of the specified Vector by a scalar.
///
/// @param v Vector to be divided.
/// @param s scalar to divide the vector by.
/// @return Vector containing the result.
/// @related Vector
template <class T, int d>
inline Vector<T, d> operator/(const Vector<T, d>& v, const T& s) {
  MATHFU_VECTOR_OPERATOR(v.data_[i] / s);
}

/// @brief Add a scalar to each element of a Vector.
///
/// @param s scalar to add to each element of a Vector.
/// @param v Vector to add the scalar to.
/// @return Vector containing the result.
/// @related Vector
template <class T, int d>
inline Vector<T, d> operator+(const T& s, const Vector<T, d>& v) {
  MATHFU_VECTOR_OPERATOR(v.data_[i] + s);
}

/// @brief Subtract a scalar from each element of a Vector.
///
/// @param s scalar to subtract from each element of a Vector.
/// @param v Vector to subtract the scalar from.
/// @return Vector containing the result.
/// @related Vector
template <class T, int d>
inline Vector<T, d> operator-(const T& s, const Vector<T, d>& v) {
  MATHFU_VECTOR_OPERATOR(v.data_[i] - s);
}

/// @brief Multiply a vector by another Vector.
///
/// In line with GLSL, this performs component-wise multiplication.
/// @param lhs First vector to use as a starting point.
/// @param rhs Second vector to multiply by.
/// @return A new Vector containing the result.
template <class T, int d>
inline Vector<T, d> operator*(const Vector<T, d>& lhs,
                              const Vector<T, d>& rhs) {
  return HadamardProductHelper(lhs, rhs);
}

/// @brief Divide a vector by another Vector.
///
/// In line with GLSL, this performs component-wise division.
/// @param lhs First vector to use as a starting point.
/// @param rhs Second vector to divide by.
/// @return A new Vector containing the result.
template <class T, int d>
inline Vector<T, d> operator/(const Vector<T, d>& lhs,
                              const Vector<T, d>& rhs) {
  MATHFU_VECTOR_OPERATOR(lhs.data_[i] / rhs[i]);
}

/// @brief Add a vector with another Vector.
///
/// @param lhs First vector to use as a starting point.
/// @param rhs Second vector to add by.
/// @return A new vector containing the result.
template <class T, int d>
inline Vector<T, d> operator+(const Vector<T, d>& lhs,
                              const Vector<T, d>& rhs) {
  MATHFU_VECTOR_OPERATOR(lhs.data_[i] + rhs[i]);
}

/// @brief subtract a vector with another Vector.
///
/// @param lhs First vector to use as a starting point.
/// @param rhs Second vector to subtract by.
/// @return A new vector containing the result.
template <class T, int d>
inline Vector<T, d> operator-(const Vector<T, d>& lhs,
                              const Vector<T, d>& rhs) {
  MATHFU_VECTOR_OPERATOR(lhs.data_[i] - rhs[i]);
}

/// @brief Multiply a vector with a scalar.
///
/// @param v Vector for the operation.
/// @param s A scalar to multiply the vector with.
/// @return A new vector containing the result.
template <class T, int d>
inline Vector<T, d> operator*(const Vector<T, d>& v, const T& s) {
  MATHFU_VECTOR_OPERATOR(v.data_[i] * s);
}

/// @brief Add a scalar to all elements of a vector.
///
/// @param v Vector for the operation.
/// @param s A scalar to add to the vector.
/// @return A new vector containing the result.
template <class T, int d>
inline Vector<T, d> operator+(const Vector<T, d>& v, const T& s) {
  MATHFU_VECTOR_OPERATOR(v.data_[i] + s);
}

/// @brief Subtract a scalar from all elements of a vector.
///
/// @param v Vector for the operation.
/// @param s A scalar to subtract from a vector.
/// @return A new vector that stores the result.
template <class T, int d>
inline Vector<T, d> operator-(const Vector<T, d>& v, const T& s) {
  MATHFU_VECTOR_OPERATOR(v.data_[i] - s);
}

/// @brief Multiply (in-place) a vector with another Vector.
///
/// In line with GLSL, this performs component-wise multiplication.
/// @param lhs First vector to use as a starting point.
/// @param rhs Second vector to multiply by.
/// @return A reference to the input <b>v</b> vector.
template <class T, int d>
inline Vector<T, d>& operator*=(Vector<T, d>& lhs, const Vector<T, d>& rhs) {
  MATHFU_VECTOR_OPERATION(lhs.data_[i] *= rhs[i]);
  return lhs;
}

/// @brief Divide (in-place) a vector by another Vector.
///
/// In line with GLSL, this performs component-wise division.
/// @param lhs First vector to use as a starting point.
/// @param rhs Second vector to divide by.
/// @return A reference to the input <b>v</b> vector.
template <class T, int d>
inline Vector<T, d>& operator/=(Vector<T, d>& lhs, const Vector<T, d>& rhs) {
  MATHFU_VECTOR_OPERATION(lhs.data_[i] /= rhs[i]);
  return lhs;
}

/// @brief Add (in-place) a vector with another Vector.
///
/// @param lhs First vector to use as a starting point.
/// @param rhs Second vector to add.
/// @return A reference to the input <b>v</b> vector.
template <class T, int d>
inline Vector<T, d>& operator+=(Vector<T, d>& lhs, const Vector<T, d>& rhs) {
  MATHFU_VECTOR_OPERATION(lhs.data_[i] += rhs[i]);
  return lhs;
}

/// @brief Subtract (in-place) another Vector from a vector.
///
/// @param lhs First vector to use as a starting point.
/// @param rhs Second vector to subtract by.
/// @return A reference to the input <b>v</b> vector.
template <class T, int d>
inline Vector<T, d>& operator-=(Vector<T, d>& lhs, const Vector<T, d>& rhs) {
  MATHFU_VECTOR_OPERATION(lhs.data_[i] -= rhs[i]);
  return lhs;
}

/// @brief Multiply (in-place) each element of a vector with a scalar.
///
/// @param v Vector for the operation.
/// @param s A scalar to multiply the vector with.
/// @return A reference to the input <b>v</b> vector.
template <class T, int d>
inline Vector<T, d>& operator*=(Vector<T, d>& v, const T& s) {
  MATHFU_VECTOR_OPERATION(v.data_[i] *= s);
  return v;
}

/// @brief Divide (in-place) each element of a vector by a scalar.
///
/// @param v Vector for the operation.
/// @param s A scalar to divide the vector by.
/// @return A reference to the input <b>v</b> vector.
template <class T, int d>
inline Vector<T, d>& operator/=(Vector<T, d>& v, const T& s) {
  MATHFU_VECTOR_OPERATION(v.data_[i] /= s);
  return v;
}

/// @brief Add (in-place) a scalar to each element of a vector.
///
/// @param v Vector for the operation.
/// @param s A scalar to add the vector to.
/// @return A reference to the input <b>v</b> vector.
template <class T, int d>
inline Vector<T, d>& operator+=(Vector<T, d>& v, const T& s) {
  MATHFU_VECTOR_OPERATION(v.data_[i] += s);
  return v;
}

/// @brief Subtract (in-place) a scalar from each element of a vector.
///
/// @param v Vector for the operation.
/// @param s A scalar to subtract from the vector.
/// @return A reference to the input <b>v</b> vector.
template <class T, int d>
inline Vector<T, d>& operator-=(Vector<T, d>& v, const T& s) {
  MATHFU_VECTOR_OPERATION(v.data_[i] -= s);
  return v;
}

/// @brief Calculate the hadamard or componentwise product of two vectors.
///
/// @param v1 First vector.
/// @param v2 Second vector.
/// @return The hadamard product of v1 and v2.
template <class T, int d>
inline Vector<T, d> HadamardProductHelper(const Vector<T, d>& v1,
                                          const Vector<T, d>& v2) {
  MATHFU_VECTOR_OPERATOR(v1[i] * v2[i]);
}

/// @brief Calculate the cross product of two vectors.
///
/// Note that this function is only defined for 3-dimensional Vectors.
/// @param v1 First vector.
/// @param v2 Second vector.
/// @return The cross product of v1 and v2.
template <class T>
inline Vector<T, 3> CrossProductHelper(const Vector<T, 3>& v1,
                                       const Vector<T, 3>& v2) {
  return Vector<T, 3>(v1[1] * v2[2] - v1[2] * v2[1],
                      v1[2] * v2[0] - v1[0] * v2[2],
                      v1[0] * v2[1] - v1[1] * v2[0]);
}

/// @brief Calculate the squared length of a vector.
///
/// @param v Vector to get the squared length of.
/// @return The length of the vector squared.
template <class T, int d>
inline T LengthSquaredHelper(const Vector<T, d>& v) {
  return DotProductHelper(v, v);
}

/// @brief Calculate the length of a vector.
///
/// @param v Vector to get the squared length of.
/// @return The length of the vector.
template <class T, int d>
inline T LengthHelper(const Vector<T, d>& v) {
  return sqrt(LengthSquaredHelper(v));
}

/// @brief Normalize a vector in-place.
///
/// @param v Vector to get the squared length of.
/// @return The length of the vector.
template <class T, int d>
inline T NormalizeHelper(Vector<T, d>& v) {
  const T length = LengthHelper(v);
  v *= (T(1) / length);
  return length;
}

/// @brief Calculate the normalized version of a vector.
///
/// @param v Vector to get the squared length of.
/// @return The normalized vector.
template <class T, int d>
inline Vector<T, d> NormalizedHelper(const Vector<T, d>& v) {
  return v * (T(1) / LengthHelper(v));
}

/// @brief Linearly interpolate two vectors.
///
/// @param v1 First vector.
/// @param v2 Second vector.
/// @param percent Percentage from v1 to v2 in range 0.0...1.0.
/// @return The hadamard product of v1 and v2.
template <class T, int d>
inline Vector<T, d> LerpHelper(const Vector<T, d>& v1, const Vector<T, d>& v2,
                               const T percent) {
  const T one_minus_percent = static_cast<T>(1.0) - percent;
  MATHFU_VECTOR_OPERATOR(one_minus_percent * v1[i] + percent * v2[i]);
}

/// @brief Generates a random vector.
///
/// The range of each component is bounded by min and max.
/// @param min Minimum value of the vector.
/// @param max Maximum value of the vector.
template <class T, int d>
inline Vector<T, d> RandomInRangeHelper(const Vector<T, d>& min,
                                        const Vector<T, d>& max) {
  Vector<T, d> result;
  MATHFU_VECTOR_OPERATION(result[i] = mathfu::RandomInRange<T>(min[i], max[i]));
  return result;
}

/// @brief Compare each component and returns max values.
///
/// @param v1 First vector.
/// @param v2 Second vector.
/// @return Max value of v1 and v2.
template <class T, int d>
inline Vector<T, d> MaxHelper(const Vector<T, d>& v1, const Vector<T, d>& v2) {
  Vector<T, d> result;
  MATHFU_VECTOR_OPERATION(result[i] = std::max(v1[i], v2[i]));
  return result;
}

/// @brief Compare each component and returns min values.
///
/// @param v1 First vector.
/// @param v2 Second vector.
/// @return Min value of v1 and v2.
template <class T, int d>
inline Vector<T, d> MinHelper(const Vector<T, d>& v1, const Vector<T, d>& v2) {
  Vector<T, d> result;
  MATHFU_VECTOR_OPERATION(result[i] = std::min(v1[i], v2[i]));
  return result;
}

/// @brief Check if val is within [range_start..range_end), denoting a
/// rectangular area.
///
/// @param val 2D vector to be tested.
/// @param range_start Starting point of the range (inclusive).
/// @param range_end Ending point of the range (non-inclusive).
/// @return Bool indicating success.
///
/// @tparam T Type of vector components to test.
template <class T>
bool InRange2D(const mathfu::Vector<T, 2>& val,
               const mathfu::Vector<T, 2>& range_start,
               const mathfu::Vector<T, 2>& range_end) {
  return InRange(val[0], range_start[0], range_end[0]) &&
         InRange(val[1], range_start[1], range_end[1]);
}

/// @cond MATHFU_INTERNAL
/// @brief Calculate the dot product of two vectors.
///
/// @param v1 First vector.
/// @param v2 Second vector.
/// @return The dot product of v1 and v2.
/// @related Vector
template <class T, int d>
static inline T DotProductHelper(const Vector<T, d>& v1,
                                 const Vector<T, d>& v2) {
  T result = 0;
  MATHFU_VECTOR_OPERATION(result += v1[i] * v2[i]);
  return result;
}
/// @endcond

/// @cond MATHFU_INTERNAL
template <class T>
static inline T DotProductHelper(const Vector<T, 2>& v1,
                                 const Vector<T, 2>& v2) {
  return v1[0] * v2[0] + v1[1] * v2[1];
}
/// @endcond

/// @cond MATHFU_INTERNAL
template <class T>
static inline T DotProductHelper(const Vector<T, 3>& v1,
                                 const Vector<T, 3>& v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}
/// @endcond

/// @cond MATHFU_INTERNAL
template <class T>
static inline T DotProductHelper(const Vector<T, 4>& v1,
                                 const Vector<T, 4>& v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3];
}
/// @endcond

template <class T, int d>
static inline bool FuzzyZeroHelper(const Vector<T, d>& v) {
    return FuzzyZero(v.LengthSquared());
}

/// @cond MATHFU_INTERNAL
template <typename T, int d, typename CompatibleT>
static inline Vector<T, d> FromTypeHelper(const CompatibleT& compatible) {
// C++11 is required for constructed unions.
#if __cplusplus >= 201103L
  // Use a union instead of reinterpret_cast to avoid aliasing bugs.
  union ConversionUnion {
    ConversionUnion() {}  // C++11.
    CompatibleT compatible;
    VectorPacked<T, d> packed;
  } u;
  static_assert(sizeof(u.compatible) == d * sizeof(T),
                "Conversion size mismatch.");

  // The read of `compatible` and write to `u.compatible` gets optimized away,
  // and this becomes essentially a safe reinterpret_cast.
  u.compatible = compatible;

  // Call the packed vector constructor with the `compatible` data.
  return Vector<T, d>(u.packed);
#else
  // Use the less-desirable memcpy technique if C++11 is not available.
  // Most compilers understand memcpy deep enough to avoid replace the function
  // call with a series of load/stores, which should then get optimized away,
  // however in the worst case the optimize away may not happen.
  // Note: Memcpy avoids aliasing bugs because it operates via unsigned char*,
  // which is allowed to alias any type.
  // See:
  // http://stackoverflow.com/questions/15745030/type-punning-with-void-without-breaking-the-strict-aliasing-rule-in-c99
  Vector<T, d> v;
  assert(sizeof(compatible) == d * sizeof(T));
  memcpy(&v, &compatible, sizeof(compatible));
  return v;
#endif  // __cplusplus >= 201103L
}
/// @endcond

/// @cond MATHFU_INTERNAL
template <typename T, int d, typename CompatibleT>
static inline CompatibleT ToTypeHelper(const Vector<T, d>& v) {
// See FromTypeHelper() for comments.
#if __cplusplus >= 201103L
  union ConversionUnion {
    ConversionUnion() {}
    CompatibleT compatible;
    VectorPacked<T, d> packed;
  } u;
  static_assert(sizeof(u.compatible) == d * sizeof(T), "Conversion size mismatch.");
  v.Pack(&u.packed);
  return u.compatible;
#else
  CompatibleT compatible;
  assert(sizeof(compatible) == d * sizeof(T));
  memcpy(&compatible, &v, sizeof(compatible));
  return compatible;
#endif  // __cplusplus >= 201103L
}
/// @endcond

/// @}

/// @addtogroup mathfu_utilities
/// @{

/// @brief Specialized version of RoundUpToPowerOf2 for vector.
template <typename T, int d>
inline Vector<T, d> RoundUpToPowerOf2(const Vector<T, d>& v) {
  Vector<T, d> ret;
  MATHFU_VECTOR_OPERATION(ret(i) = RoundUpToPowerOf2(v(i)));
  return ret;
}
/// @}

}  // namespace mathfu

// Include the specializations to avoid template errors.
// For example, if you include vector.h, use Vector<float, 3>, and then
// include vector_3.h, you the compiler will generate an error since you're
// specializing something that has already been instantiated.
#include "mathfu/internal/vector_2_simd.h"
#include "mathfu/internal/vector_3_simd.h"
#include "mathfu/internal/vector_4_simd.h"

#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif

#endif  // MATHFU_VECTOR_H_
