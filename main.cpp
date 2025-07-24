#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <iostream>
#include <fstream>

// 添加VTK相关头文件
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef CGAL::Surface_mesh<K::Point_3>                          Mesh;
typedef boost::graph_traits<Mesh>::edge_descriptor              edge_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

// 根据旋转轴和角度构造旋转变换矩阵
CGAL::Aff_transformation_3<K> create_rotation_matrix(const Vector_3& axis, double angle_rad) {
    const double x = axis.x();
    const double y = axis.y();
    const double z = axis.z();
    const double cos_a = cos(angle_rad);
    const double sin_a = sin(angle_rad);
    const double one_minus_cos_a = 1 - cos_a;

    return CGAL::Aff_transformation_3<K>(
        cos_a + x * x * one_minus_cos_a, x * y * one_minus_cos_a - z * sin_a, x * z * one_minus_cos_a + y * sin_a, 0,
        y * x * one_minus_cos_a + z * sin_a, cos_a + y * y * one_minus_cos_a, y * z * one_minus_cos_a - x * sin_a, 0,
        z * x * one_minus_cos_a - y * sin_a, z * y * one_minus_cos_a + x * sin_a, cos_a + z * z * one_minus_cos_a, 0
        );
}
//生成圆柱体
Mesh create_cylinder(const CylinderParams& params, const Point_3& base, const Vector_3& direction) {
    Mesh mesh;
    std::vector<Mesh::Vertex_index> bottom, top;

    Mesh::Vertex_index bottom_center = mesh.add_vertex(Point_3(0, 0, 0));
    Mesh::Vertex_index top_center = mesh.add_vertex(Point_3(0, 0, params.height));

    for (int i = 0; i < params.segments; ++i) {
        double angle = 2 * CGAL_PI * i / params.segments;
        double x = params.radius * cos(angle);
        double y = params.radius * sin(angle);
        bottom.push_back(mesh.add_vertex(Point_3(x, y, 0)));
        top.push_back(mesh.add_vertex(Point_3(x, y, params.height)));
    }

    for (int i = 0; i < params.segments; ++i) {
        int j = (i + 1) % params.segments;
        mesh.add_face(bottom[i], bottom[j], top[j]);
        mesh.add_face(top[j], top[i], bottom[i]);
    }

    for (int i = 0; i < params.segments; ++i) {
        int j = (i + 1) % params.segments;
        mesh.add_face(bottom_center, bottom[j], bottom[i]);
        mesh.add_face(top_center, top[i], top[j]);
    }

    // 计算旋转参数
    const Vector_3 z_axis(0, 0, 1);
    Vector_3 rot_axis = CGAL::cross_product(z_axis, direction);
    const double axis_length = CGAL::sqrt(rot_axis.squared_length());

    if (axis_length > 1e-6) {
        // 规范化旋转轴
        rot_axis = rot_axis / axis_length;
        const double cos_theta = z_axis * direction / CGAL::sqrt(direction.squared_length());
        const double angle = acos(cos_theta);

        // 创建旋转矩阵
        auto rot = create_rotation_matrix(rot_axis, angle);
        PMP::transform(rot, mesh);
    }
    CGAL::Aff_transformation_3<K> trans(CGAL::TRANSLATION, base - CGAL::ORIGIN); PMP::transform(trans, mesh);

    // 缝合圆柱边界
    PMP::stitch_borders(mesh);
    return mesh;
}
//生成环形分布点
std::vector<Point_3> generate_base_points(const Mesh& mesh, const DrainageConfig& config, Point_3& center_out) {
    auto bbox = PMP::bbox(mesh);
    Point_3 center(
        (bbox.xmin() + bbox.xmax()) / 2,
        (bbox.ymin() + bbox.ymax()) / 2,
        bbox.zmin() + config.base_height
    );
    center_out = center;

    double radius = config.radius_scale * std::min(bbox.xmax() - bbox.xmin(), bbox.ymax() - bbox.ymin()) / 2;
    std::vector<Point_3> positions;
    for (int i = 0; i < config.num_holes; ++i) {
        double angle = 2 * CGAL_PI * i / config.num_holes;
        positions.emplace_back(
            center.x() + radius * cos(angle),
            center.y() + radius * sin(angle),
            center.z()
        );
    }
    return positions;
}


