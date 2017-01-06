####################################
#Imports
from libcpp.string cimport string
from libcpp cimport bool as bool_cpp
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr

####################################
#Soplex Rational
cdef extern from "rational.h" namespace "soplex":
    cdef cppclass Rational:
        Rational() except +
        #double operator double() not yet supported
        Rational(const double& r) except +


####################################
#IPO Vector


cdef extern from "ipo/vectors.h" namespace "ipo":
    cdef cppclass MutableVector:
        MutableVector() except +
    cdef cppclass Vector:

        Vector() except +
        void swap(Vector& other)

        #Vector& operator=(const Vector& other)
        #Vector& operator=(const MutableVector& other)
        #not possible in python to overload


    cdef cppclass ReferenceCountedVector:
        size_t size()
        size_t index(size_t position)
        bool_cpp isSorted()
        double approximation(size_t position)
        const Rational& value(size_t position)
        bool_cpp operator==(const ReferenceCountedVector& other)
        bool_cpp operator!=(const ReferenceCountedVector& other)
        bool_cpp operator<(const ReferenceCountedVector& other)

    cdef struct InnerDescription:
        vector[Vector] points;
        vector[Vector] rays;

####################################
#IPO Linear Constraint

cdef extern from "ipo/linear_constraint.h" namespace "ipo":
    cdef cppclass LinearConstraint:
        LinearConstraint() except +
        LinearConstraint(char type, const Vector& normal, const Rational& rhs) except +

        bool_cpp operator==(const LinearConstraint& other)
        bool_cpp operator<(const LinearConstraint& other)
        bool_cpp isEquation()
        char type()
        const Vector& normal()
        const Rational& rhs()
        Rational getMaximumNorm()

        bool_cpp definesCompleteFace()
        bool_cpp definesEmptyFace()
        bool_cpp definesTrivialFace()

        int evaluatePoint(const Vector& point)
        int evaluateRay(const Vector& ray)

    ctypedef vector[LinearConstraint] AffineOuterDescription

####################################
#IPO Space

cdef extern from "<iostream>" namespace "std":
    cdef cppclass ostream:
        ostream()
    ostream cout

cdef extern from "ipo/space.h" namespace "ipo":
    cdef cppclass Space:
        Space() except +
        Space(const Space& other) except +

        bool_cpp operator==(const Space& other)
        bool_cpp operator!=(const Space& other)
        const string& operator[](size_t variable)

        size_t dimension()

        void printVector(ostream& stream, const Vector& vector)
        void printLinearForm(ostream& stream, const Vector& linearForm)
        void printLinearConstraint(ostream& stream, const LinearConstraint& constraint)

####################################
#IPO Polyhedron

cdef extern from "ipo/oracles.h" namespace "ipo":
    ctypedef size_t HeuristicLevel
    cdef cppclass OracleBase:
        pass

cdef extern from "ipo/polyhedron.h" namespace "ipo":
    cdef cppclass Polyhedron:

        cppclass Face:
            Face(const LinearConstraint& inequality) except +

            const LinearConstraint& inequality()
            bool_cpp hasDimension()
            const AffineOuterDescription& outerDescription()
            const InnerDescription& innerDescription()
            int dimension()

        Polyhedron(const shared_ptr[OracleBase]& oracle) except +

        Space space()
        size_t numPoints()
        size_t numRays()
        size_t numInequalities()
        shared_ptr[Face] inequalityToFace(const LinearConstraint& constraint)
        #void affineHull(shared_ptr[Face]& face, vector[AffineHullHandler*]& handlers)
        void affineHull(shared_ptr[Face]& face)
        #void affineHull(vector[AffineHullHandler*]& handlers)
        void affineHull()
        int dimension()
        const AffineOuterDescription& affineHullOuterDescription()
        const InnerDescription& affineHullInnerDescription()
        void setAffineHullLastCheapHeuristicLevel(HeuristicLevel lastCheapHeuristicLevel)
        void setAffineHullLastModerateHeuristicLevel(HeuristicLevel lastModerateHeuristicLevel)
        void addConstraint(const LinearConstraint& constraint)
        void getFaces(vector[shared_ptr[Face]]& constraints, bool_cpp onlyInequalities, bool_cpp onlyWithDimension)


####################################
#IPO ScipOracle Python Wrapper

cdef extern from "ipo/python_wrapper.h" namespace "ipo":
    cdef cppclass OracleControllerBase:
        OracleControllerBase()

    cdef cppclass ScipOracleController(OracleControllerBase):
        ScipOracleController(string name) except +
        ScipOracleController(string name, OracleControllerBase* next) except +
        ScipOracleController(string name, OracleControllerBase* next, const shared_ptr[MixedIntegerSet] mixedIntegerSet) except +
        #ScipOracleController(string name, const shared_ptr[MixedIntegerSet] mixedIntegerSet) except +


        string name()
        int heuristicLevel_ScipOracle()
        int heuristicLevel_CacheOracle()
        shared_ptr[OracleBase] getCacheOracle()

        InnerDescription affineHullInner(int outputMode)
        AffineOuterDescription affineHullOuter(int outputMode)

####################################
#IPO MixedIntegerSet class

cdef extern from "ipo/mip.h" namespace "ipo":


    cdef cppclass MixedIntegerSet:
        #constructor needs scip instance
        #MixedIntegerSet(SCIP* scip)
        const Space& space()
        size_t numVariables()
        size_t numRows()
        const Variable& variable(size_t variable)
        const LinearConstraint& rowConstraint(size_t row)
        const vector[LinearConstraint]& rowConstraints()
        const string& rowName(size_t row)
        bool_cpp isIntegral(size_t variable)
        void setFace(const LinearConstraint& newFace)
cdef extern from "ipo/mip.h" namespace "ipo::MixedIntegerSet":
    cdef struct Variable:
        bool_cpp integral
        Rational upperBound
        Rational lowerBound

cdef extern from "ipo/python_wrapper.h" namespace "ipo":
    const shared_ptr[MixedIntegerSet] getMixedIntegerSet(string filename)

####################################
#IPO Test

cdef extern from "ipo/test.h" namespace "ipo":
    cdef cppclass Foo:
        Foo() except +
        string printFoo()
