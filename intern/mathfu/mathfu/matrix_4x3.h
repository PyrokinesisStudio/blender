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
#ifndef MATHFU_MATRIX_4X3_H_
#define MATHFU_MATRIX_4X3_H_

#include "mathfu/matrix.h"

/// @file mathfu/matrix_4x3.h
/// @brief Matrix class and functions.
/// @addtogroup mathfu_matrix
///
/// MathFu provides a generic Matrix implementation which is specialized
/// for 4x4 matrices to take advantage of optimization opportunities using
/// SIMD instructions.


namespace mathfu {

template <class T>
class Matrix<T, 4, 3> {
 public:
  /// @brief Construct a Matrix of uninitialized values.
  inline Matrix<T, 4, 3>() {}

  /// @brief Construct a Matrix from another Matrix copying each element.
  ////
  /// @param m Matrix that the data will be copied from.
  inline Matrix<T, 4, 3>(const Matrix<T, 4, 3>& m) {
    MATHFU_MAT_OPERATION(data_[i] = m.data_[i]);
  }

  /// @brief Construct a Matrix from a single float.
  ///
  /// @param s Scalar value used to initialize each element of the matrix.
  explicit inline Matrix<T, 4, 3>(const T& s) {
    MATHFU_MAT_OPERATION((data_[i] = Vector<T, 4>(s)));
  }

  /// @brief Creates a Matrix from twelve floats.
  ///
  /// @note This method only works with Matrix<float, 4, 3>.
  ///
  ///
  /// @param s00 Value of the first row and column.
  /// @param s10 Value of the second row, first column.
  /// @param s20 Value of the third row, first column.
  /// @param s30 Value of the fourth row, first column.
  /// @param s01 Value of the first row, second column.
  /// @param s11 Value of the second row and column.
  /// @param s21 Value of the third row, second column.
  /// @param s31 Value of the fourth row, second column.
  /// @param s02 Value of the first row, third column.
  /// @param s12 Value of the second row, third column.
  /// @param s22 Value of the third row and column.
  /// @param s32 Value of the fourth row, third column.
  inline Matrix<T, 4, 3>(const T& s00, const T& s10, const T& s20, const T& s30,
						 const T& s01, const T& s11, const T& s21, const T& s31,
						 const T& s02, const T& s12, const T& s22, const T& s32) {
    data_[0] = Vector<T, 4>(s00, s10, s20, s30);
    data_[1] = Vector<T, 4>(s01, s11, s21, s31);
    data_[2] = Vector<T, 4>(s02, s12, s22, s32);
  }

  /// @brief Create 4x4 Matrix from 4, 4 element vectors.
  ///
  /// @note This method only works with a 4x4 Matrix.
  ///
  /// @param column0 Vector used for the first column.
  /// @param column1 Vector used for the second column.
  /// @param column2 Vector used for the third column.
  /// @param column3 Vector used for the fourth column.
  inline Matrix<T, 4, 3>(const Vector<T, 4>& column0, const Vector<T, 4>& column1,
						 const Vector<T, 4>& column2) {
    data_[0] = column0;
    data_[1] = column1;
    data_[2] = column2;
  }

  /// @brief Create a Matrix from the first row * column elements of an array.
  ///
  /// @param a Array of values that the matrix will be iniitlized to.
  explicit inline Matrix<T, 4, 3>(const T* const a) {
    MATHFU_MAT_OPERATION((data_[i] = Vector<T, 4>(&a[i * 3])));
  }

  /// @brief Create a Matrix from the first row * column elements of an array.
  ///
  /// @param a Array of values that the matrix will be iniitlized to.
  explicit inline Matrix<T, 4, 3>(const T a[3][4]) {
    MATHFU_MAT_OPERATION((data_[i] = Vector<T, 4>(a[i])));
  }

  /// @brief Create a Matrix from an array of "3", "4" element packed
  /// vectors.
  ///
  /// @param vectors Array of "3", "4" element packed vectors.
  explicit inline Matrix<T, 4, 3>(const VectorPacked<T, 4>* const vectors) {
    MATHFU_MAT_OPERATION((data_[i] = Vector<T, 4>(vectors[i])));
  }

  /// @brief Access an element of the matrix.
  ///
  /// @param row Index of the row to access.
  /// @param column Index of the column to access.
  /// @return Const reference to the element.
  inline const T& operator()(const int row, const int column) const {
    return data_[column][row];
  }

  /// @brief Access an element of the Matrix.
  ///
  /// @param row Index of the row to access.
  /// @param column Index of the column to access.
  /// @return Reference to the data that can be modified by the caller.
  inline T& operator()(const int row, const int column) {
    return data_[column][row];
  }

  /// @brief Access an element of the Matrix.
  ///
  /// @param i Index of the element to access in flattened memory.  Where
  /// the column accessed is i / 4 and the row is i % 4.
  /// @return Reference to the data that can be modified by the caller.
  inline const T& operator()(const int i) const { return operator[](i); }

  /// @brief Access an element of the Matrix.
  ///
  /// @param i Index of the element to access in flattened memory.  Where
  /// the column accessed is i / 4 and the row is i % 4.
  /// @return Reference to the data that can be modified by the caller.
  inline T& operator()(const int i) { return operator[](i); }

  /// @brief Access an element of the Matrix.
  ///
  /// @param i Index of the element to access in flattened memory.  Where
  /// the column accessed is i / 4 and the row is i % 4.
  /// @return Const reference to the data.
  inline const T& operator[](const int i) const {
    return const_cast<Matrix<T, 4, 3>*>(this)->operator[](i);
  }

  /// @brief Access an element of the Matrix.
  ///
  /// @param i Index of the element to access in flattened memory.  Where
  /// the column accessed is i / 4 and the row is i % 4.
  /// @return Reference to the data that can be modified by the caller.
  inline T& operator[](const int i) {
#if defined(MATHFU_COMPILE_WITH_PADDING)
    // In this case Vector<T, 3> is padded, so the element offset must be
    // accessed using the array operator.
    if (4 == 3) {
      const int row = i % 4;
      const int col = i / 4;
      return data_[col][row];
    } else {
      return reinterpret_cast<T*>(data_)[i];
    }
#else
    return reinterpret_cast<T*>(data_)[i];
#endif  // defined(MATHFU_COMPILE_WITH_PADDING)
  }

  /// @brief Pack the matrix to an array of "4" element vectors,
  /// one vector per matrix column.
  ///
  /// @param vector Array of "3" entries to write to.
  inline void Pack(VectorPacked<T, 4>* const vector) const {
    MATHFU_MAT_OPERATION(GetColumn(i).Pack(&vector[i]));
  }

  /// @brief Pack the matrix to an array of "4" element vectors,
  /// one vector per matrix column.
  ///
  /// @param vector Array of "3" entries to write to.
  inline void Pack(T a[3][4]) const {
    MATHFU_MAT_OPERATION(GetColumn(i).Pack(a[i]));
  }

  inline void Pack(T a[3 * 4]) const {
    MATHFU_MAT_OPERATION(GetColumn(i).Pack(&a[i * 3]));
  }

  /// @cond MATHFU_INTERNAL
  /// @brief Access a column vector of the Matrix.
  ///
  /// @param i Index of the column to access.
  /// @return Reference to the data that can be modified by the caller.
  inline Vector<T, 4>& GetColumn(const int i) { return data_[i]; }

  /// @brief Access a column vector of the Matrix.
  ///
  /// @param i Index of the column to access.
  /// @return Const reference to the data.
  inline const Vector<T, 4>& GetColumn(const int i) const {
    return data_[i];
  }
  /// @endcond

  /// @brief Negate this Matrix.
  ///
  /// @return Matrix containing the result.
  inline Matrix<T, 4, 3> operator-() const {
    MATHFU_MAT_OPERATOR(-data_[i]);
  }

  /// @brief Add a Matrix to this Matrix.
  ///
  /// @param m Matrix to add to this Matrix.
  /// @return Matrix containing the result.
  inline Matrix<T, 4, 3> operator+(
      const Matrix<T, 4, 3>& m) const {
    MATHFU_MAT_OPERATOR(data_[i] + m.data_[i]);
  }

  /// @brief Subtract a Matrix from this Matrix.
  ///
  /// @param m Matrix to subtract from this Matrix.
  /// @return Matrix containing the result.
  inline Matrix<T, 4, 3> operator-(
      const Matrix<T, 4, 3>& m) const {
    MATHFU_MAT_OPERATOR(data_[i] - m.data_[i]);
  }

  /// @brief Add a scalar to each element of this Matrix.
  ///
  /// @param s Scalar to add to this Matrix.
  /// @return Matrix containing the result.
  inline Matrix<T, 4, 3> operator+(const T& s) const {
    MATHFU_MAT_OPERATOR(data_[i] + s);
  }

  /// @brief Subtract a scalar from each element of this Matrix.
  ///
  /// @param s Scalar to subtract from this matrix.
  /// @return Matrix containing the result.
  inline Matrix<T, 4, 3> operator-(const T& s) const {
    MATHFU_MAT_OPERATOR(data_[i] - s);
  }

  /// @brief Multiply each element of this Matrix with a scalar.
  ///
  /// @param s Scalar to multiply with this Matrix.
  /// @return Matrix containing the result.
  inline Matrix<T, 4, 3> operator*(const T& s) const {
    MATHFU_MAT_OPERATOR(data_[i] * s);
  }

  /// @brief Divide each element of this Matrix with a scalar.
  ///
  /// @param s Scalar to divide this Matrix with.
  /// @return Matrix containing the result.
  inline Matrix<T, 4, 3> operator/(const T& s) const {
    return (*this) * (1 / s);
  }

  /// @brief Multiply this Matrix with another Matrix.
  ///
  /// @param m Matrix to multiply with this Matrix.
  /// @return Matrix containing the result.
  inline Matrix<T, 4, 3> operator*(
      const Matrix<T, 4, 3>& m) const {
    Matrix<T, 4, 3> result;
    TimesHelper(*this, m, &result);
    return result;
  }

  /// @brief Add a Matrix to this Matrix (in-place).
  ///
  /// @param m Matrix to add to this Matrix.
  /// @return Reference to this class.
  inline Matrix<T, 4, 3>& operator+=(
      const Matrix<T, 4, 3>& m) {
    MATHFU_MAT_SELF_OPERATOR(data_[i] += m.data_[i]);
  }

  /// @brief Subtract a Matrix from this Matrix (in-place).
  ///
  /// @param m Matrix to subtract from this Matrix.
  /// @return Reference to this class.
  inline Matrix<T, 4, 3>& operator-=(
      const Matrix<T, 4, 3>& m) {
    MATHFU_MAT_SELF_OPERATOR(data_[i] -= m.data_[i]);
  }

  /// @brief Add a scalar to each element of this Matrix (in-place).
  ///
  /// @param s Scalar to add to each element of this Matrix.
  /// @return Reference to this class.
  inline Matrix<T, 4, 3>& operator+=(const T& s) {
    MATHFU_MAT_SELF_OPERATOR(data_[i] += s);
  }

  /// @brief Subtract a scalar from each element of this Matrix (in-place).
  ///
  /// @param s Scalar to subtract from each element of this Matrix.
  /// @return Reference to this class.
  inline Matrix<T, 4, 3>& operator-=(const T& s) {
    MATHFU_MAT_SELF_OPERATOR(data_[i] -= s);
  }

  /// @brief Multiply each element of this Matrix with a scalar (in-place).
  ///
  /// @param s Scalar to multiply with each element of this Matrix.
  /// @return Reference to this class.
  inline Matrix<T, 4, 3>& operator*=(const T& s) {
    MATHFU_MAT_SELF_OPERATOR(data_[i] *= s);
  }

  /// @brief Divide each element of this Matrix by a scalar (in-place).
  ///
  /// @param s Scalar to divide this Matrix by.
  /// @return Reference to this class.
  inline Matrix<T, 4, 3>& operator/=(const T& s) {
    return (*this) *= (1 / s);
  }

  /// @brief Multiply this Matrix with another Matrix (in-place).
  ///
  /// @param m Matrix to multiply with this Matrix.
  /// @return Reference to this class.
  inline Matrix<T, 4, 3>& operator*=(
      const Matrix<T, 4, 3>& m) {
    const Matrix<T, 4, 3> copy_of_this(*this);
    TimesHelper(copy_of_this, m, this);
    return *this;
  }

  inline Vector<T, 4> GetEuler() const {
     return EulerHelper(*this);
  }

  /// @brief Calculate the inverse of this Matrix.
  ///
  /// This calculates the inverse Matrix such that
  /// <code>(m * m).Inverse()</code> is the identity.
  /// @return Matrix containing the result.
  inline Matrix<T, 4, 3> Inverse() const {
    Matrix<T, 4, 3> inverse;
    InverseHelper<false>(*this, &inverse);
    return inverse;
  }

  /// @brief Calculate the inverse of this Matrix.
  ///
  /// This calculates the inverse Matrix such that
  /// <code>(m * m).Inverse()</code> is the identity.
  /// By contrast to Inverse() this returns whether the matrix is invertible.
  ///
  /// The invertible check simply compares the calculated determinant with
  /// Constants<T>::GetDeterminantThreshold() to roughly determine whether the
  /// matrix is invertible.  This simple check works in common cases but will
  /// fail for corner cases where the matrix is a combination of huge and tiny
  /// values that can't be accurately represented by the floating point
  /// datatype T.  More extensive checks (relative to the input values) are
  /// possible but <b>far</b> more expensive, complicated and difficult to
  /// test.
  /// @return Whether the matrix is invertible.
  inline bool InverseWithDeterminantCheck(
      Matrix<T, 4, 3>* const inverse) const {
    return InverseHelper<true>(*this, inverse);
  }

  /// @brief Calculate the transpose of this Matrix.
  ///
  /// @return The transpose of the specified Matrix.
  inline Matrix<T, 3, 4> Transpose() const {
    Matrix<T, 3, 4> transpose;
    MATHFU_UNROLLED_LOOP(
        i, 3, MATHFU_UNROLLED_LOOP(
                        j, 4, transpose.GetColumn(j)[i] = GetColumn(i)[j]))
    return transpose;
  }

  /// @brief Get the 2-dimensional translation of a 2-dimensional affine
  /// transform.
  ///
  /// @note 2-dimensional affine transforms are represented by 3x3 matrices.
  /// @return Vector with the first two components of column 2 of this Matrix.
  inline Vector<T, 2> TranslationVector2D() const {
    MATHFU_STATIC_ASSERT(4 == 3 && 3 == 3);
    return Vector<T, 2>(data_[2][0], data_[2][1]);
  }

  /// @brief Get the 3-dimensional translation of a 3-dimensional affine
  /// transform.
  ///
  /// @note 3-dimensional affine transforms are represented by 4x4 matrices.
  /// @return Vector with the first three components of column 3.
  inline Vector<T, 3> TranslationVector3D() const {
    MATHFU_STATIC_ASSERT(4 == 4 && 3 == 4);
    return Vector<T, 3>(data_[3][0], data_[3][1], data_[3][2]);
  }

  /// @brief Load from any byte-wise compatible external matrix.
  ///
  /// Format should be `3` vectors, each holding `4` values of type T.
  ///
  /// Use this for safe conversion from external matrix classes.
  /// Often, external libraries will have their own matrix types that are,
  /// byte-for-byte, exactly the same as mathfu::Matrix. This function allows
  /// you to load a mathfu::Matrix from those external types, without potential
  /// aliasing bugs that are caused by casting.
  ///
  /// @note If your external type gives you access to a T*, then you can
  ///       equivalently use the Matrix(const T*) constructor.
  ///
  /// @param compatible reference to a byte-wise compatible matrix structure;
  ///                   array of 3 x 4 Ts.
  /// @returns `compatible` loaded as a mathfu::Matrix.
  template <typename CompatibleT>
  static inline Matrix<T, 4, 3> FromType(const CompatibleT& compatible) {
    return FromTypeHelper<T, 4, 3, CompatibleT>(compatible);
  }

  /// @brief Load into any byte-wise compatible external matrix.
  ///
  /// Format should be `3` vectors, each holding `4` values of type T.
  ///
  /// Use this for safe conversion to external matrix classes.
  /// Often, external libraries will have their own matrix types that are,
  /// byte-for-byte, exactly the same as mathfu::Matrix. This function allows
  /// you to load an external type from a mathfu::Matrix, without potential
  /// aliasing bugs that are caused by casting.
  ///
  /// @param m reference to mathfu::Matrix to convert.
  /// @returns CompatibleT loaded from m.
  template <typename CompatibleT>
  static inline CompatibleT ToType(const Matrix<T, 4, 3>& m) {
    return ToTypeHelper<T, 4, 3, CompatibleT>(m);
  }

  /// @brief Calculate the outer product of two Vectors.
  ///
  /// @return Matrix containing the result.
  static inline Matrix<T, 4, 3> OuterProduct(
      const Vector<T, 4>& v1, const Vector<T, 3>& v2) {
    return OuterProductHelper(v1, v2);
  }

  /// @brief Calculate the hadamard / component-wise product of two matrices.
  ///
  /// @param m1 First Matrix.
  /// @param m2 Second Matrix.
  /// @return Matrix containing the result.
  static inline Matrix<T, 4, 3> HadamardProduct(
      const Matrix<T, 4, 3>& m1, const Matrix<T, 4, 3>& m2) {
    MATHFU_MAT_OPERATOR(m1[i] * m2[i]);
  }

  /// @brief Calculate the identity Matrix.
  ///
  /// @return Matrix containing the result.
  static inline Matrix<T, 4, 3> Identity() {
    return IdentityHelper<T, 4, 3>();
  }

  /// @brief Multiply a Vector by a Matrix.
  ///
  /// @param v Vector to multiply.
  /// @param m Matrix to multiply.
  /// @return Matrix containing the result.
  friend inline Vector<T, 3> operator*(
      const Vector<T, 4>& v, const Matrix<T, 4, 3>& m) {
    const int d = 3;
    MATHFU_VECTOR_OPERATOR((Vector<T, 4>::DotProduct(m.data_[i], v)));
  }

  // Dimensions of the matrix.
  /// Number of 4 in the matrix.
  static const int kRows = 4;
  /// Number of 3 in the matrix.
  static const int kColumns = 3;
  /// Total number of elements in the matrix.
  static const int kElements = 4 * 3;

  MATHFU_DEFINE_CLASS_SIMD_AWARE_NEW_DELETE

 private:
  Vector<T, 4> data_[3];
};
/// @}

#endif  // MATHFU_MATRIX_4X3_H_

