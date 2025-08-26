import plotly.graph_objects as go


def visualize_retinal(df, double_bonds, carbon_annotations, cutoff_distance=1.6):
    fig = go.Figure()

    for idx, row in df.iterrows():
        fig.add_trace(go.Scatter3d(
            x=[row['x']],
            y=[row['y']],
            z=[row['z']],
            mode='markers+text',
            marker=dict(size=5, color='blue'),
            text=[carbon_annotations.get(idx, '') if idx in carbon_annotations.keys() else 'CH3'],
            textposition='bottom center'
        ))

    def are_neighbors(index1, index2):
        x1, y1, z1 = df.iloc[index1]
        x2, y2, z2 = df.iloc[index2]
        distance = ((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) ** 0.5
        return distance <= cutoff_distance

    # Add single and double bonds
    for i in range(len(df)):
        for j in range(i + 1, len(df)):
            if are_neighbors(i, j):
                if (i, j) in double_bonds or (j, i) in double_bonds:
                    # Add a double bond (thicker line)
                    fig.add_trace(go.Scatter3d(
                        x=[df.iloc[i]['x'], df.iloc[j]['x']],
                        y=[df.iloc[i]['y'], df.iloc[j]['y']],
                        z=[df.iloc[i]['z'], df.iloc[j]['z']],
                        mode='lines',
                        line=dict(color='red', width=5),
                        name='Double Bond'
                    ))
                else:
                    # Add a single bond (thinner line)
                    fig.add_trace(go.Scatter3d(
                        x=[df.iloc[i]['x'], df.iloc[j]['x']],
                        y=[df.iloc[i]['y'], df.iloc[j]['y']],
                        z=[df.iloc[i]['z'], df.iloc[j]['z']],
                        mode='lines',
                        line=dict(color='black', width=2),
                        name='Single Bond'
                    ))

    # Update layout for a better view
    fig.update_layout(
        scene=dict(
            xaxis_title='X Axis',
            yaxis_title='Y Axis',
            zaxis_title='Z Axis'
        ),
        showlegend=False
    )
    fig.show()