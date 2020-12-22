#pragma once

#include <memory>
#include <vector>


/*! @class Field
 *
 *  @brief A multidimensional array representing a field
 *
 *	This class is an abstraction of a field.
 *	Its number of dimensions and dimension sizes are completely free.
 *	It can represent both steady and unsteady fields.
 *	The field represented can be scalar, vectorial, matricial or tensorial.
 *	The field can be of any spacial dimension too.
 *	This class allows to manage easily dynamically allocated, contiguous, arbitrarily big and high-dimensional arrays.
 *	The user is responsible for the definition and meaning of each dimension.
 */
class Field
{
private:
	size_t m_nb_dim;
	std::vector<size_t> m_dimensions;

	std::unique_ptr<double[]> ptr_data;

public:
	/*! @brief The field constructor
	 *	
	 *	This constructor takes the size of each dimension as a std::vector<size_t> and dynamically allocates the required memory.
	 *	It creates a Field object which number of dimensions is equal to the length of the std::vector.
	 *
	 *  @param dimensions A std::vector<size_t> indicating the size of each dimension.
	 *	It is often useful to implicitely declare this std::vector with an initialization list.
	 *
	 *	@warning The constructor allocates by does not initialize the array.
	 */
	Field(std::vector<size_t> dimensions = {});

	/*! @brief The field destructor
	 *
	 *	This destructor first deallocates the memory used by the instance, and then destroys the instance.
	 */
	~Field();

	/*! @brief A method reshaping the instance
	 *	
	 *	This method takes the size of each dimension as a std::vector<size_t> and dynamically allocates the required memory.
	 *	It reshapes the Field object to match the number of dimensions to the length of the std::vector.
	 *
	 *  @param dimensions A std::vector<size_t> indicating the size of each dimension.
	 *	It is often useful to implicitely declare this std::vector with an initialization list.
	 *
	 *	@warning This method allocates by does not initialize the array.
	 *
	 *	@warning Any previously stored data will be lost by calling this method.
	 */
	void reshape(std::vector<size_t> dimensions);
	
	/*! @brief The field getter/setter
	 *
	 *	Depending on the context, this operator acts both as a getter and a setter for the base elements of the field (double).
	 *	Its behavior is similar to the [] operator on classic arrays, but is multidimentional.
	 *
	 *  @param index A std::vector<size_t> which length equals the number of dimensions of the field, and indicating which element to access.
	 *	It is often useful to implicitely declare this std::vector with an initialization list.
	 *
	 *  @returns A reference to the chosen element of the field.
	 *
	 *  @warning This operator does not check for out of bounds requests, which would result in undefined behavior.
	 *
	 *  @bug Due to the absence of boundary checks, this operator should not be used on the degenerate case of a field of dimension zero. 
	 */
	double& operator()(std::vector<size_t> index);

	/*! @brief A getter to the memory address of the first element of the field
	 *
	 *	@returns A pointer to the first element of the array
	 *
	 *	@warning This method is meant to be used only by the Python buffer protocol and should not be called by the user directly.
	 */
	double* get_ptr() const;

	/*! @brief A getter to the memory size of each element of the field
	 *
	 *	@returns The size of an element of the array
	 *
	 *	@warning This method is meant to be used only by the Python buffer protocol and should not be called by the user directly.
	 */
	size_t get_element_size() const;

	/*! @brief A getter to the number of dimensions of the field
	 *
	 *	@returns The number of dimensions of the field
	 */
	size_t get_nb_dimensions() const;

	/*! @brief A getter to the size of each dimension of the field
	 *
	 *	@returns A std::vector containing the size of each dimension
	 */
	std::vector<size_t> get_dimensions_sizes() const;

	/*! @brief A getter to the adjusted sizes of the elements of the field
	 *
	 *	@returns The apparent size of an element for each dimension
	 *
	 *	@warning This method is meant to be used only by the Python buffer protocol and should not be called by the user directly.
	 */
	std::vector<size_t> get_adjusted_element_sizes() const;
};