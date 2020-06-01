#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtCharts>

QT_BEGIN_NAMESPACE
using namespace QtCharts;

namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = nullptr);
	~MainWindow();
	void connect_axes(int a, int b);
	QScatterSeries* create_ser_without_resistance( double V_x, double V_y, double g);
	QLineSeries* create_numer_approx ( int i, const int N, double _h, double k, double m, double g, double k1, double k2,
							  double k3, double k4, double l1, double l2, double l3, double l4, double pV_x[],
							  double pV_y[], double pY[], double pV[], double pX[]);
	QLineSeries* create_relativism( double k1, double l1, double k2, double l2, double k3, double l3, double k4,
						   double l4, double _h, double pV_y[], double pV_x[], double pX[], double pY[],
						   double pV[], double m, double k, double g, int j, const int N);
	QLineSeries* create_department_mass ( double x0, double y0, double m_p, double k, double dt, double V_xp,
										  double V_yp, double m_p2, double g, double V_t);


private:
	Ui::MainWindow *ui;

	QChart *chrt;
	QChart *chart;

	QValueAxis *axisX, *axisY;
};
#endif // MAINWINDOW_H
