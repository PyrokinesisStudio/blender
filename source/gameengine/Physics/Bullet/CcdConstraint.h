#ifndef __CCD_CONSTRAINT_H__
#define __CCD_CONSTRAINT_H__

#include "PHY_IConstraint.h"

class btTypedConstraint;

class CcdConstraint : public PHY_Constraint
{
private:
	btTypedConstraint *m_constraint;

public:
	CcdConstraint(btTypedConstraint *constraint, int id, PHY_ConstraintType type);
	virtual ~CcdConstraint();

	virtual void SetParam(int param, float value, float value1);
	virtual float GetParam(int param);

	virtual int GetIdentifier() const;
	virtual PHY_ConstraintType GetType() const;
};

#endif  // __CCD_CONSTRAINT_H__
