#define int PetesInt
#define float PetscScalr

class FEI {
public:
	int n, row, *column;
	float *value, rightItem;
	FEI(int n, int row) : n(n), row(row) {
		column = new int[n];
		value = new float[n];
	}
	void set(int column, float value) {
		this->column[column] = value;
	}
}