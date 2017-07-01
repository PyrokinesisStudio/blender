#ifndef MATHFU_FRUSTUM_H_
#define MATHFU_FRUSTUM_H_

#include "mathfu/matrix.h"
#include "mathfu/vector.h"

namespace mathfu {

static const Vector<float, 3> normalizedBox[8] = {
	Vector<float, 3>(-1.0f, -1.0f, -1.0f),
	Vector<float, 3>(-1.0f, 1.0f, -1.0f),
	Vector<float, 3>(1.0f, 1.0f, -1.0f),
	Vector<float, 3>(1.0f, -1.0f, -1.0f),
	Vector<float, 3>(-1.0f, -1.0f, 1.0f),
	Vector<float, 3>(-1.0f, 1.0f, 1.0f),
	Vector<float, 3>(1.0f, 1.0f, 1.0f),
	Vector<float, 3>(1.0f, -1.0f, 1.0f)
};

inline void FrustumBox(const Matrix<float, 4>& m, std::array<Vector<float, 3>, 8>& box)
{
	for (unsigned short i = 0; i < 8; ++i) {
		const Vector<float, 3>& p3 = normalizedBox[i];
		const Vector<float, 4> p4 = m * Vector<float, 4>(p3, 1.0f);

		box[i] = p4.xyz() / p4.w;
	}
}

inline void FrustumAabb(const Matrix<float, 4>& m, Vector<float, 3>& min, Vector<float, 3>& max)
{
	for (unsigned short i = 0; i < 8; ++i) {
		const Vector<float, 3>& p3 = normalizedBox[i];
		const Vector<float, 4> p4 = m * Vector<float, 4>(p3, 1.0f);
		const Vector<float, 3> co = p4.xyz() / p4.w;

		if (i == 0) {
			min = max = co;
		}
		else {
			min = Vector<float, 3>::Min(min, co);
			max = Vector<float, 3>::Max(max, co);
		}
	}
}

}

#endif  // MATHFU_FRUSTUM_H_
