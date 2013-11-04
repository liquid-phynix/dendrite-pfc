#include <set>
#include <list>
#include <vector>
#include <utility>

struct Edge {
  Edge(int _row, int _col, bool _ES):
    row(_row), col(_col), index(0), ES(_ES){}
  int row, col;
  mutable int index;
  bool ES; // true = east, false = south
};

struct EdgeCompare {
  bool operator() (const Edge& lhs, const Edge& rhs) const {
    if(lhs.row == rhs.row){
      if(lhs.col == rhs.col){
        return lhs.ES < rhs.ES;
      } else return lhs.col < rhs.col;
    } else return lhs.row < rhs.row;
  }
};

struct Point {
  double x, y;
  Point(): x(0), y(0){}
  //  Point(Edge& e){ y = e.row; x = e.col; }
  Point(double _x, double _y){ y = _y; x = _x; }
  Point add(Point other){ return Point(x + other.x, y + other.y); }
  Point smul(double s){ return Point(s * x, s * y); }
  friend std::ostream& operator<<(std::ostream& out, Point& obj){ out << obj.x << "\t" << obj.y; return out; }
  static Point make_point(double cont, int r1, int c1, double v1, int r2, int c2, double v2){
    Point zero(c1, r1); double t;
    t = (v1 - cont) / (v1 - v2);
    return zero.smul(1 - t).add(Point(c2, r2).smul(t));
  }
};

void contour(Dist_Matrix<double>& mat, double cont, const char* fmt, int it){
  const Edge case1_1(0, 0, 0),  case1_2(1, 0, 1);
  const Edge case2_1(1, 0, 1),  case2_2(0, 1, 0);
  const Edge case3_1(0, 0, 0),  case3_2(0, 1, 0);
  const Edge case4_1(0, 0, 1),  case4_2(0, 1, 0);
  const Edge case5_1(0, 0, 0),  case5_2(0, 0, 1),  case5_3(1, 0, 1), case5_4(0, 1, 0);
  const Edge case6_1(0, 0, 1),  case6_2(1, 0, 1);
  const Edge case7_1(0, 0, 0),  case7_2(0, 0, 1);
  const Edge case8_1(0, 0, 0),  case8_2(0, 0, 1);
  const Edge case9_1(0, 0, 1),  case9_2(1, 0, 1);
  const Edge case10_1(0, 0, 0), case10_2(1, 0, 1), case10_3(0, 0, 1), case10_4(0, 1, 0);
  const Edge case11_1(0, 0, 1), case11_2(0, 1, 0);
  const Edge case12_1(0, 0, 0), case12_2(0, 1, 0);
  const Edge case13_1(1, 0, 1), case13_2(0, 1, 0);
  const Edge case14_1(0, 0, 0), case14_2(1, 0, 1);
  
  std::set<Edge, EdgeCompare> edges;
  std::list<std::pair<std::set<Edge, EdgeCompare>::iterator, std::set<Edge, EdgeCompare>::iterator> > points;
  std::pair<std::set<Edge, EdgeCompare>::iterator, bool> ret1, ret2;

  int count = 0;

  const int rows = mat.local_rows() - 1, cols = mat.local_cols() - 1;
  for(int r = 0; r < rows; r++){
    for(int c = 0; c < cols; c++){
      bool D = mat.at(r,   c)   > cont;

      count += D;

      bool C = mat.at(r,   c+1) > cont;
      bool A = mat.at(r+1, c)   > cont;
      bool B = mat.at(r+1, c+1) > cont;
      switch(D * 8 + C * 4 + B * 2 + A){
      case 0: break;
      case 1:
        ret1 = edges.insert(Edge(r + case1_1.row, c + case1_1.col, case1_1.ES));
        ret2 = edges.insert(Edge(r + case1_2.row, c + case1_2.col, case1_2.ES));
        points.push_back(std::make_pair(ret1.first, ret2.first));
        break;
      case 2:
        ret1 = edges.insert(Edge(r + case2_1.row, c + case2_1.col, case2_1.ES));
        ret2 = edges.insert(Edge(r + case2_2.row, c + case2_2.col, case2_2.ES));
        points.push_back(std::make_pair(ret1.first, ret2.first));
        break;
      case 3:
        ret1 = edges.insert(Edge(r + case3_1.row, c + case3_1.col, case3_1.ES));
        ret2 = edges.insert(Edge(r + case3_2.row, c + case3_2.col, case3_2.ES));
        points.push_back(std::make_pair(ret1.first, ret2.first));
        break;
      case 4:
        ret1 = edges.insert(Edge(r + case4_1.row, c + case4_1.col, case4_1.ES));
        ret2 = edges.insert(Edge(r + case4_2.row, c + case4_2.col, case4_2.ES));
        points.push_back(std::make_pair(ret1.first, ret2.first));
        break;
      case 5:
        ret1 = edges.insert(Edge(r + case5_1.row, c + case5_1.col, case5_1.ES));
        ret2 = edges.insert(Edge(r + case5_2.row, c + case5_2.col, case5_2.ES));
        points.push_back(std::make_pair(ret1.first, ret2.first));
        ret1 = edges.insert(Edge(r + case5_3.row, c + case5_3.col, case5_3.ES));
        ret2 = edges.insert(Edge(r + case5_4.row, c + case5_4.col, case5_4.ES));
        points.push_back(std::make_pair(ret1.first, ret2.first));
        break;
      case 6:
        ret1 = edges.insert(Edge(r + case6_1.row, c + case6_1.col, case6_1.ES));
        ret2 = edges.insert(Edge(r + case6_2.row, c + case6_2.col, case6_2.ES));
        points.push_back(std::make_pair(ret1.first, ret2.first));
        break;
      case 7:
        ret1 = edges.insert(Edge(r + case7_1.row, c + case7_1.col, case7_1.ES));
        ret2 = edges.insert(Edge(r + case7_2.row, c + case7_2.col, case7_2.ES));
        points.push_back(std::make_pair(ret1.first, ret2.first));
        break;
      case 8:
        ret1 = edges.insert(Edge(r + case8_1.row, c + case8_1.col, case8_1.ES));
        ret2 = edges.insert(Edge(r + case8_2.row, c + case8_2.col, case8_2.ES));
        points.push_back(std::make_pair(ret1.first, ret2.first));
        break;
      case 9:
        ret1 = edges.insert(Edge(r + case9_1.row, c + case9_1.col, case9_1.ES));
        ret2 = edges.insert(Edge(r + case9_2.row, c + case9_2.col, case9_2.ES));
        points.push_back(std::make_pair(ret1.first, ret2.first));
        break;
      case 10:
        ret1 = edges.insert(Edge(r + case10_1.row, c + case10_1.col, case10_1.ES));
        ret2 = edges.insert(Edge(r + case10_2.row, c + case10_2.col, case10_2.ES));
        points.push_back(std::make_pair(ret1.first, ret2.first));
        ret1 = edges.insert(Edge(r + case10_3.row, c + case10_3.col, case10_3.ES));
        ret2 = edges.insert(Edge(r + case10_4.row, c + case10_4.col, case10_4.ES));
        points.push_back(std::make_pair(ret1.first, ret2.first));
        break;
      case 11:
        ret1 = edges.insert(Edge(r + case11_1.row, c + case11_1.col, case11_1.ES));
        ret2 = edges.insert(Edge(r + case11_2.row, c + case11_2.col, case11_2.ES));
        points.push_back(std::make_pair(ret1.first, ret2.first));
        break;
      case 12:
        ret1 = edges.insert(Edge(r + case12_1.row, c + case12_1.col, case12_1.ES));
        ret2 = edges.insert(Edge(r + case12_2.row, c + case12_2.col, case12_2.ES));
        points.push_back(std::make_pair(ret1.first, ret2.first));
        break;
      case 13:
        ret1 = edges.insert(Edge(r + case13_1.row, c + case13_1.col, case13_1.ES));
        ret2 = edges.insert(Edge(r + case13_2.row, c + case13_2.col, case13_2.ES));
        points.push_back(std::make_pair(ret1.first, ret2.first));
        break;
      case 14:
        ret1 = edges.insert(Edge(r + case14_1.row, c + case14_1.col, case14_1.ES));
        ret2 = edges.insert(Edge(r + case14_2.row, c + case14_2.col, case14_2.ES));
        points.push_back(std::make_pair(ret1.first, ret2.first));
        break;
      case 15: break;
      default: break; }}}

  std::list<Edge> edges_list;
  int idx = 0;
  for(auto it = edges.begin(); it != edges.end(); it++){
    it->index = idx;
    edges_list.push_back(*it);
    idx++; }
  // for(auto it = points.begin(); it != points.end(); it++){
  //   std::cout << "p1: " << it->first->index << ", p2: " << it->second->index << std::endl;
  // }

  int sendbuf[2] = { (int) points.size(), (int) edges_list.size() };
  int recvbuf[2];
  MPI::COMM_WORLD.Scan(sendbuf, recvbuf, 2, MPI::INT, MPI::SUM);
  int points_offset = recvbuf[0] - points.size();
  int edges_offset = recvbuf[1] - edges_list.size();

//  std::cout << "rank " << MPI::COMM_WORLD.Get_rank() << " reports " << edges_list.size() << " points and offset " << edges_offset << " and count " << count << std::endl;

  std::vector<Point> real_points; real_points.resize(edges_list.size());
  std::vector<std::pair<int, int> > line_endpoints; line_endpoints.resize(points.size());

  { int idx = 0;
    for(auto it = points.begin(); it != points.end(); it++, idx++){
      line_endpoints[idx].first = it->first->index + edges_offset;
      line_endpoints[idx].second = it->second->index + edges_offset; }}

  
  { int idx = 0;
    int offset = mat.get_rows_start();
    for(auto it = edges_list.begin(); it != edges_list.end(); it++, idx++){
      int dr = it->ES ? 0 : 1;
      int dc = it->ES ? 1 : 0;
      real_points[idx] =
        Point::make_point(cont,
                          it->row + offset,      it->col,      mat.at(it->row, it->col),
                          it->row + offset + dr, it->col + dc, mat.at(it->row + dr, it->col + dc)); }}

  char fn[128];
  { // saving points
    sprintf(fn, fmt, "points", it);
    LAST { 
      std::fstream fs(fn,  std::fstream::out | std::fstream::trunc);
      fs << "# rows=" << recvbuf[1] << " cols=2"
         << " dtype=" << NUMPY_ENDIAN << "f8" << std::endl; }
    MPI::COMM_WORLD.Barrier();
    MPI::File fh = MPI::File::Open(MPI::COMM_WORLD, fn, MPI::MODE_WRONLY | MPI::MODE_APPEND, MPI::INFO_NULL);
    fh.Set_errhandler(MPI::ERRORS_ARE_FATAL);
    fh.Set_view(fh.Get_position() + edges_offset * 2 * sizeof(double), MPI::DOUBLE, MPI::DOUBLE, "native", MPI::INFO_NULL);
    fh.Write_all(real_points.data(), 2 * real_points.size(), MPI::DOUBLE);
    fh.Close();
    LAST std::cout << "points saved to '" << fn << "'" << std::endl; }

  { // saving line segment endpoint indices
    sprintf(fn, fmt, "endpoint_indices", it);
    LAST { 
      std::fstream fs(fn,  std::fstream::out | std::fstream::trunc);
      fs << "# rows=" << recvbuf[0] << " cols=2"
         << " dtype=" << NUMPY_ENDIAN << "i4" << std::endl; }
    MPI::COMM_WORLD.Barrier();
    MPI::File fh = MPI::File::Open(MPI::COMM_WORLD, fn, MPI::MODE_WRONLY | MPI::MODE_APPEND, MPI::INFO_NULL);
    fh.Set_errhandler(MPI::ERRORS_ARE_FATAL);
    fh.Set_view(fh.Get_position() + points_offset * 2 * sizeof(int), MPI::INT, MPI::INT, "native", MPI::INFO_NULL);
    fh.Write_all(line_endpoints.data(), 2 * line_endpoints.size(), MPI::INT);
    fh.Close();
    LAST std::cout << "endpoints saved to '" << fn << "'" << std::endl; }
}
