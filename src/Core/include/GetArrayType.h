#ifndef _GETArrayTYPE_H
#define _GETArrayTYPE_H
#define UNKNOWN_TYPE 999
namespace Core
{
	template <typename TYPE, unsigned int M>
	class ArrayDefault;
	//This class gets the integer type of a matrix. The default implementation returns type unknown
	template<typename Array>
	class GetArrayType 
	{
	public:
		inline explicit GetArrayType() { }
		static unsigned int ArrayType() { return UNKNOWN_TYPE; }
	};

	template<>
	class GetArrayType<ArrayDefault<double,9999> >
	{
	public:
		inline explicit GetArrayType() { }
		static unsigned int ArrayType() { return 1; }
	};

	template<>
	class GetArrayType<ArrayDefault<float,9999> >
	{
	public:
		inline explicit GetArrayType() { }
		static unsigned int ArrayType() { return 2; }
	};

	template<>
	class GetArrayType<ArrayDefault<double,6> >
	{
	public:
		inline explicit GetArrayType() { }
		static unsigned int ArrayType() { return 4; }
	};
}
#endif
